import numpy as np
import xarray as xr
import xarray_regrid


eau_precip=0
list_precip=[]



def mass(ds: xr.Dataset, var_name: str = "concentration") -> float:
    """Masse physique = sum(C * dz), avec dz = z_top - z_bot."""
    dz = (ds["z_top"] - ds["z_bot"]).to_numpy()
    C = ds[var_name].to_numpy()
    return float(np.sum(C * dz))


z_top_ref = np.array([2000,4000.0,6000, 8000.0,10000,12000.0])  # 3 mailles de 4 km
ds_ref = make_vertical_ds(z_top_ref, var_name="concentration", init="top_cell_100")
dz = 2000

print(ds_ref)

data_0 = ds_ref
z_top_last = data_0['z_top'].values[-1]


def advect_down(ds: xr.Dataset, V: float, dt: float) -> xr.Dataset:
    """
    Décale les mailles vers le bas à vitesse V (m/s) pendant dt (s).
    On décale:
    - la coord z (centres)
    - z_top et z_bot
    """
    shift = -V * dt
    out = ds.copy(deep=True)
    out = out.assign_coords(level=out["level"] + shift)
    out["z_top"] = out["z_top"] + shift
    out["z_bot"] = out["z_bot"] + shift
    return out




data_bin_i_t = data_0
list_precip.append(data_0["concentration"].values)
# On itère sur le nombre de bin

dt = 1
list_data = []
for k in range(6):
    for diam in range(1):

        speed = 3000

        data_i_dt_tempo = advect_down(data_bin_i_t, V=speed, dt=dt)

        nb_stitche_create = 1

        z_mid_last_i = data_i_dt_tempo["level"].values[-1]

        new_grid = np.sort(np.unique(np.append(data_i_dt_tempo["level"].values,z_mid_last_i+dz)))

        data_i_dt_tempo = data_i_dt_tempo.reindex(level = new_grid)

        

        data_i_dt_tempo["z_top"] = data_i_dt_tempo["z_top"].copy(data=data_i_dt_tempo["z_top"].data).fillna(z_mid_last_i+3/2*dz)
        data_i_dt_tempo["z_bot"] = data_i_dt_tempo["z_bot"].copy(data=data_i_dt_tempo["z_bot"].data).fillna(z_mid_last_i+1/2*dz)
        data_i_dt_tempo["concentration"] = data_i_dt_tempo["concentration"].copy(data=data_i_dt_tempo["concentration"].data).fillna(0)

        nb_stitche_create  = 1

        nb_stitche_create  += 1


        while dz + data_i_dt_tempo["level"].values[-1] <=z_top_ref[-1]:

            

            z_mid_last_i = data_i_dt_tempo["level"].values[-1]

            new_grid = np.sort(np.unique(np.append(data_i_dt_tempo["level"].values,z_mid_last_i+dz)))

            data_i_dt_tempo = data_i_dt_tempo.reindex(level = new_grid)

            

            data_i_dt_tempo["z_top"] = data_i_dt_tempo["z_top"].copy(data=data_i_dt_tempo["z_top"].data).fillna(z_mid_last_i+3/2*dz)
            data_i_dt_tempo["z_bot"] = data_i_dt_tempo["z_bot"].copy(data=data_i_dt_tempo["z_bot"].data).fillna(z_mid_last_i+1/2*dz)
            data_i_dt_tempo["concentration"] = data_i_dt_tempo["concentration"].copy(data=data_i_dt_tempo["concentration"].data).fillna(0)

            nb_stitche_create  += 1






        ds_on_old = data_i_dt_tempo.regrid.conservative(data_bin_i_t, time_dim=None).compute()

        ds_on_old["z_top"] = data_0["z_top"].copy(data=data_0["z_top"].data)
        ds_on_old["z_bot"] = data_0["z_bot"].copy(data=data_0["z_bot"].data)
        ds_on_old["concentration"] = ds_on_old["concentration"].copy(data=ds_on_old["concentration"].data.todense())


        M0 = mass(data_bin_i_t)
        M1 = mass(ds_on_old)
        
        eau_precip += M0-M1



        data_bin_i_t = ds_on_old


        list_data.append(ds_on_old["concentration"].values)
        list_precip.append(eau_precip)

print(list_data,list_precip)