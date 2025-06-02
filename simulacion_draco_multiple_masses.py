import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import quantity_support
from gammapy.maps import MapAxis, WcsGeom
from gammapy.irf import load_irf_dict_from_file
from gammapy.datasets import MapDataset, Datasets
from gammapy.data import Observation
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.modeling.models import (
    SkyModel, PointSpatialModel, FoVBackgroundModel, Models
)
from gammapy.astro.darkmatter import DarkMatterAnnihilationSpectralModel
from gammapy.astro.darkmatter.profiles import EinastoProfile
from gammapy.modeling import select_nested_models
from gammapy.stats.utils import ts_to_sigma
from gammapy.estimators import ParameterEstimator
from concurrent.futures import ProcessPoolExecutor, as_completed
import time


def run_simulation(i, mass_dm, dataset_base, channel, jfact_draco, center):
    spatial_model = PointSpatialModel(
        lon_0=center.ra,
        lat_0=center.dec,
        frame="icrs"
    )

    spatial_model.lon_0.min = spatial_model.lon_0.value - 0.5
    spatial_model.lon_0.max = spatial_model.lon_0.value + 0.5
    spatial_model.lat_0.min = spatial_model.lat_0.value - 0.5
    spatial_model.lat_0.max = spatial_model.lat_0.value + 0.5

    spectral_model = DarkMatterAnnihilationSpectralModel(
        mass=mass_dm, channel=channel, jfactor=jfact_draco
    )

    spectral_model.scale.scan_min = 1
    spectral_model.scale.scan_max = 1e4
    spectral_model.scale.scan_n_values = 10
    spectral_model.scale.frozen = False

    dataset = dataset_base.copy()
    sky_model = SkyModel(spatial_model=spatial_model, spectral_model=spectral_model, name="draco-dm")
    bkg_model = FoVBackgroundModel(dataset_name=dataset.name)
    dataset.models = Models([sky_model, bkg_model])
    dataset.fake()

    if i == 0:
        quantity_support()
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        dataset.counts.sum_over_axes().plot(ax=axes[0], add_cbar=True)
        axes[0].set_title("Cuentas simuladas")
        dataset.npred_background().sum_over_axes().plot(ax=axes[1], add_cbar=True)
        axes[1].set_title("Fondo esperado")
        dataset.npred_signal().sum_over_axes().plot(ax=axes[2], add_cbar=True)
        axes[2].set_title("Modelo DM")
        plt.tight_layout()
        plt.savefig(f"mapas_mass_{int(mass_dm.value)}TeV_rep_{i}.png")
        plt.close()

    datasets = Datasets([dataset])
    datasets.models = dataset.models
    parameters = [datasets.models["draco-dm"].parameters["scale"]]
    null_values = [0]
    n_free = len(sky_model.parameters.free_parameters)

    nested = select_nested_models(
        datasets=datasets,
        parameters=parameters,
        null_values=null_values,
        n_sigma=-np.inf,
        n_free_parameters=n_free
    )
    significance = ts_to_sigma(nested["ts"], df=n_free)

    estimator = ParameterEstimator(
        n_sigma=2,
        n_sigma_ul=2,
        selection_optional="all",
        reoptimize=True
    )
    result = estimator.run(datasets, parameter="scale")
    ul = result.get("scale_ul", np.nan)

    return {
        "mass_TeV": mass_dm.value,
        "repetition": i,
        "significance": significance,
        "scale_UL": ul,
        "UL_<sigmav>_95CL": ul * 3e-26 if ul is not None else np.nan
    }


def main():
    parser = argparse.ArgumentParser(description="Simulación de Draco para una masa de DM.")
    parser.add_argument("--mass", type=float, required=True, help="Masa de DM en TeV")
    parser.add_argument("--n_repeat", type=int, default=10, help="Número de repeticiones por masa")
    parser.add_argument("--channel", type=str, default="h", help="Canal de aniquilación (default=h)")
    args = parser.parse_args()

    mass_dm = args.mass * u.TeV
    n_repeat = args.n_repeat
    channel = args.channel

    jfact_draco = 1.66e18 * u.GeV**2 / u.cm**5
    draco_pos = SkyCoord(260.05, 57.915, frame="icrs", unit="deg")
    rho_s_msun_kpc3 = 1.3e7 * (u.M_sun / u.kpc**3)
    rho_s_GeV_cm3 = rho_s_msun_kpc3.to(u.GeV / u.cm**3, equivalencies=u.mass_energy())
    draco_profile = EinastoProfile(r_s=0.91 * u.kpc, rho_s=rho_s_GeV_cm3)

    energy_edges = np.logspace(-1, 1, 10)
    energy_reco = MapAxis.from_edges(energy_edges, unit="TeV", name="energy", interp="log")
    energy_true = MapAxis.from_edges(energy_edges, unit="TeV", name="energy_true", interp="log")

    geom_draco = WcsGeom.create(
        skydir=draco_pos,
        binsz=0.05,
        width=1.0,
        frame="icrs",
        axes=[energy_reco]
    )

    irfs_path = r"C:\Users\pauli\GAMMAPY-DATA\cta-prod5-zenodo-fitsonly-v0.1\fits\Prod5-North-20deg-NorthAz-4LSTs09MSTs.180000s-v0.1.fits"
    irfs = load_irf_dict_from_file(irfs_path)

    empty = MapDataset.create(
        geom=geom_draco,
        name="dataset-simu-draco",
        energy_axis_true=energy_true
    )

    pointing = draco_pos
    livetime = 500 * u.h
    time_ref = Time("2025-06-17T01:00:00", format="isot", scale="utc")
    obs = Observation.create(pointing=pointing, livetime=livetime, irfs=irfs, reference_time=time_ref)

    maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
    mask_maker = SafeMaskMaker(methods=["offset-max"], offset_max=2.5 * u.deg)

    dataset_base = maker.run(empty.copy(), obs)
    dataset_base = mask_maker.run(dataset_base, obs)

    center = dataset_base.geoms["geom"].center_skydir

    results = []

    print(f"\nLanzando {n_repeat} simulaciones para masa = {mass_dm.value:.1f} TeV...\n")
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [
            executor.submit(run_simulation, i, mass_dm, dataset_base, channel, jfact_draco, center)
            for i in range(n_repeat)
        ]
        for idx, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)

            elapsed = time.time() - start_time
            avg_time = elapsed / idx
            remaining = avg_time * (n_repeat - idx)

            print(f"[{idx}/{n_repeat}] Simulación repetición #{result['repetition']} completada.")
            print(f"   Tiempo transcurrido: {elapsed:.1f}s | Tiempo estimado restante: {remaining:.1f}s")

    df = pd.DataFrame(results)

    print("\nResultados de la simulación:")
    print(df)

    outname = f"results_mass_{int(mass_dm.value)}TeV_{channel}.csv"
    df.to_csv(outname, index=False)
    print(f"\nResultados guardados en: {outname}")

    total_time = time.time() - start_time
    print(f"\n Tiempo total de simulación: {total_time:.1f} segundos")


if __name__ == "__main__":
    main()

