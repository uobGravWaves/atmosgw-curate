# Supplying data

What files the loaders read, how to lay them out, and where each product
comes from. Output formats are in `docs/design/data_formats.md`.

## Where the loaders look

Each dataset has a folder (`datasets.toml`; override it in your settings
file, see the README, or pass `file_path=` to a loader). Some datasets have a
second folder for planetary-wave-filtered files (`pws`, used with
`hindley_pws=True`), and AIRS has one for its 2-D product (`2d`); keep those
apart from the instrument's own files.

Below a dataset's folder, any layout works: a flat folder of files (e.g.
`~/HIRDLS/*.he5`), year folders, or any other nesting. For each day (or
month, year or granule) a loader looks for the file's name:

- the date-named readers (HIRDLS, MLS, SOFIE, MIPAS daily; SABER monthly) search the whole
  folder for a name containing the date, as their MATLAB originals do;
- the others (ACE, GNSS, planetary-wave files, AIRS, IASI) try the usual path
  first (e.g. `yyyy/ace_v52_yyyydddd.mat`, `yyyy/ddd/airs_yyyy_ddd_ggg.nc`),
  then search for that name anywhere below the folder (one folder listing per
  load: under a second for the largest folders checked, AIRS and IASI, about
  50 000 entries each);
- the newer readers (LIMS, ROMEX, Aeolus, DWD sondes, IAGOS, FAAM) search
  the whole folder for their file names; IAGOS also looks inside zip
  archives there;
- ALOMAR files are looked for in the folder itself; Sampler MAT-files in
  `sample_path` itself, its netCDF files anywhere below it.

If more than one file has the name or date looked for, the first in search
order is used (at each level, the contents of sub-folders before the files
beside them, each in name order) and the
choice is warned and listed in the output attribute `multiple_matches`; the
newest MLS version is chosen instead of the first (and the DWD, IAGOS and
FAAM readers choose as described below). A file at its usual path is
used without looking further, except for IASI, whose files of a day are
gathered from its day folder and anywhere else below the IASI folder (files
in another day's folder are left out). Do not keep two copies of a product (e.g. two
versions) under one folder unless you want this choice made for you.

File names must keep the patterns below: the date is taken from the name.

## The products

For each product: the file names (dates `yyyy` year, `ddd` day of year,
`mm` month, `dd` day of month), the file type, the variables read, and where
the files come from. "Checked" is the version the readers were tested on;
other versions may work but have not been checked. The version each file
states is recorded in `ds.attrs["source_versions"]`. Fill values are as in the
files checked; every reader also masks the fill value or `missing_value` a
netCDF file declares.

### Limb sounders (`load_limb_sounders`)

**HIRDLS** (Aura)

- Files: one per day, any name containing `_yyyydddd`, HDF-EOS5
  (e.g. `HIRDLS-Aura_L2_v07-00-20-c01_2007d001.he5`).
- Variables: swath `HIRDLS`, `Geolocation Fields/` `Latitude`, `Longitude`,
  `Time` (s since 1993-01-01), `Altitude` (m), `Pressure` (hPa),
  `ScanUpFlag`; `Data Fields/` `Temperature` (K), `TemperaturePrecision`
  (negative: a priori dominated); [profiles x levels]; fill -999. The other
  fields of the product can be read with `additional_vars`
  (`_limb_hirdls.GEOLOCATION_FIELDS`, `DATA_FIELDS`).
- Checked: HIRDLS Level 2 V7 (PGEVersion V1.37.8).
- Provenance: NASA GES DISC, DOI 10.5067/AURA/HIRDLS/DATA201 (stated in the
  files).

**MLS** (Aura)

- Files: one per day, any name containing `_yyyydddd`, HDF-EOS5 Level 2 GP
  temperature (e.g. `MLS-Aura_L2GP-Temperature_v06-01-c01_2007d001.he5`).
- Variables: swath `Temperature`, `L2gpValue` (K), `L2gpPrecision`,
  `Quality`, `Status`, `Convergence`; geolocation `Latitude`, `Longitude`,
  `Time` (s since 1993-01-01), `Pressure` (hPa, one set of levels); the
  `_FillValue` each field declares.
- Checked: v6.01 (and v5 files, where both are present the newest is used).
  The screening rules follow each version's data quality document
  (`_limb_mls.screen`).
- Provenance: NASA GES DISC; v6 DOI 10.5067/AURA/MLS/DATA2620 (stated in the
  files; v5 files state their own).

**SABER** (TIMED)

- Files: one per month, any name containing `<Month><yyyy>` (any case),
  netCDF (e.g. `SABER_Temp_O3_H2O_January2020_v2.0.nc`).
- Variables: [events x levels] `ktemp` (K), `pressure` (mbar), `tpaltitude`
  (km), `tplatitude`, `tplongitude`, `time` (ms since midnight of `date`);
  per event `date` (yyyyddd); `missing_value` -999.
- Checked: Level 2A Temp_O3_H2O files v2.0, Data_Product_Version 02.07 and
  02.08.
- Provenance: GATS Inc. for the SABER team (https://saber.gats-inc.com);
  monthly custom files.

**SOFIE** (AIM)

- Files: one per day, any name containing `Level2_yyyyddd`, netCDF
  (e.g. `SOFIE_Level2_2007152_01.3.nc`).
- Variables: [events x altitudes] `Temperature` (K), `Pressure` (mbar);
  `Altitude` (km, one set); per event `Latitude_83km`, `Longitude_83km`,
  `Time_UT` (hours of the day; the files say "msec"); `missing_value` -1e24.
- Checked: Level 2 V1.3 (DP_Version 01.3).
- Provenance: GATS Inc. for the SOFIE team (https://sofie.gats-inc.com).

**MIPAS** (Envisat)

- Files: one per day, any name containing `_yyyymmdd_nom`, netCDF.
- Variables: [profiles x levels] `TEM` (K), `PRE` (hPa), `hgt` (km); per
  profile `latitude`, `longitude`, `time` (hours since the start of the
  day); `_FillValue` -999.
- Checked: Oxford retrieval v1.40, nominal mode.
- Provenance: University of Oxford MIPAS retrievals, distributed by CEDA
  (catalogue record 0ce969a5aef0fb1a7b00d7136706b76c).

**ACE-FTS** (SCISAT)

- Files: one per day, `ace_v52_yyyydddd.mat`.
- Variables: [profiles x levels] `Lat`, `Lon`, `Time` (MATLAB datenums),
  `Alt` (km), `Temp` (K), `Pres` (atmospheres: converted to hPa); NaN for
  missing.
- Checked: v5.2 retrievals (from the file name).
- Provenance: placeholder, to be updated later: these MAT-files are a
  conversion of the ACE-FTS Level 2 v5.2 data; who made them, from which
  files, and how to obtain or regenerate them is not yet documented.

**GNSS radio occultation**

- Files: one per day, `merged_ro_yyyy_ddd.mat`.
- Variables: [profiles x levels] `Temp` (K), `Pres` (Pa: converted to hPa),
  `Lat`, `Lon` (NaN rows where a retrieval gives no per-level positions),
  `PTRes` (original vertical resolution, units not stated); `MSL_alt` (km,
  one row for all profiles); a table `MetaData` with per-profile `occid`
  (ending `yyyymmddHHMM`), `latitude`, `longitude`, `time` (datenum),
  `AtmSource` (the source file name, giving centre and retrieval version);
  optional `Settings.Version` (merge version). Fill 9.97e36.
- Checked: files merging JPL, UCAR and ROM SAF retrievals, 2001-2025.
- Provenance: placeholder, to be updated later: these files merge several
  processing centres' retrievals; who made them, from which products, and how
  to obtain or regenerate them is not yet documented.

**ALOMAR lidar** (Andoya, Norway)

- Files: one per year, `yyyy_V2026-07-27_file4juli.h5`, HDF5, directly in
  the ALOMAR folder.
- Variables: one group per measurement session, each with `T1_lidar` (K),
  `T1_lidar_unc` (K) [times x heights], `time_lidar` (POSIX s), `z_lidar`
  (km).
- Checked: version V2026-07-27.
- Provenance: placeholder, to be updated later.

**Planetary-wave-filtered files** ("Hindley23" format, `hindley_pws=True`)

- Files: one per day, `<INSTRUMENT>_PWs_yyyydddd.mat` (e.g.
  `SABER_PWs_2007d001.mat`) in the instrument's `pws` folder.
- Variables: [profiles x levels] `Lat`, `Lon`, `Time` (datenums), `Alt`
  (km), `Temp_Residual` (K), optional `Pres` (hPa; else estimated from
  height); `Temp_PW` (K) [profiles x levels x modes]; optional `Note` (code
  version).
- Provenance: placeholder, to be updated later.

**Sampler output** (`Sampler`)

- Files: `sampled_<datenum>.mat` in `sample_path` (e.g.
  `sampled_733043.mat` for 2007-01-01).
- Variables: struct `Sampled_Data` with [levels x profiles] `Alt`, `Lat`,
  `Lon`, `Time`, `T`, `Prs`, `SourceProf`.
- Provenance: made by running the Wright and Hindley (2018) sampler,
  https://github.com/corwin365/Model-Sampling.

Where a day has no `sampled_<datenum>.mat`, the first netCDF file below
`sample_path` whose name ends `yyyydddd.nc` is read: the sampler's netCDF
output, as in the SWANS DYAMOND runs (e.g.
`SWANS/HIRDLS/dyamond_geos3km/HIRDLS-dyamond_geos3km-2020d020.nc`):
[levels x profiles] `Lat`, `Lon`, `Time` (datenums), `Alt`, `Prs`, `Temp`,
`ObsProf` (and `TSimple`, `Observed`, `ObsFile`, not read); fill -9999
(`_FillValue`). Checked: SWANS HIRDLS, MLS, GNSS and SABER (GEOS 3 km),
2020-01-20. Pass `sampled_from` to name the model in the output (else "a
model (not named)").

**LIMS** (Nimbus-7)

- Files: one per day, `LIMS_2_yyyymmdd_V6.ASCII` (V6 revised Level 2,
  ASCII). Each file holds complete orbits, so it runs into the next day
  (and two start before their own day: 1979-03-27, 1979-05-15); the
  previous and next day's files are also read.
- Variables: per scan a fixed-width header line (read by column: the
  fields can run together), with the 30 km tangent point `alat`, `alon`
  (degrees east), `iday` (day of the year), the time (GMT; columns read as
  `hour`, `minute`, `second`), `idn`, `iud`,
  `iad`; then 6 channel lines; then 109 level lines of 18 values (`zt_std`
  km, `pt_std` mb, `tt_std` K, radiances, gradients, mixing ratios,
  `geo_heights` km); fill 1e24. `zt_std` is geometric height and
  `geo_heights` geopotential (hypsometric check, 1978-11-08).
- Checked: the 201 files of 1978-10-25 to 1979-05-28 on eepc-0184 (every
  header in range; values of the User's Guide example reproduced; every
  one of the 630,053 scans returned once over the 214 days they cover). The
  files of 1978-12-01 to 12-06 also have `idn` 0, which the guide does not
  define (1 day, 2 night). The cloud and removed-scan lists named in the
  README are not held.
- Provenance: LIMS V6 data (GATS, NASA GES DISC); User's Guide and README
  beside the data.

**ROMEX** (GNSS radio occultation)

- Files: `romex_<mission>rmx[...]_100m_v1a[...].nc`, one or more per
  mission, 2022-09-01 to 2022-11-30; netCDF.
- Variables: `t` (s since 1970-01-01; 0: no time), `altitude` (500 levels,
  10-59.9 km), [profile x level] `lat_full`, `lon_full`, `temp` (degrees C),
  `press` (mb), `ref`, `Bend_ang`; per profile `lat`, `lon`, `azim`,
  `is_bad`, `filestamp`. No global attributes.
- `altitude` is geometric height (not stated in the files; hypsometric
  and ERA5 checks, `docs/validation.md`).
- Checked: the 15 files on eepc-0184, every day 2022-09-01 to 11-30
  (3,126,854 profiles, about 34,000 a day); against the merged GNSS
  files (same occultations, median differences -0.26 to +0.08 K) and ERA5.
  About 0.1% of profiles not marked `is_bad` exceed 330 K above 30 km
  (up to about 370 K); they are left as distributed.
- Provenance: placeholder, to be updated later: who processed these files
  and how to obtain them (ROMEX is taken to be the IROWG Radio Occultation
  Modeling Experiment: UNVERIFIED).

**Aeolus** (ALADIN Level 2B HLOS winds)

- Files: `Aeolus_L2B_monthly_yyyy_mm.nc` (made from weekly VirES downloads
  by `get_aeolus.py`, beside the data; a "monthly" file runs past its
  month, so each day is also looked for in the previous month's file).
- Variables: per observation, for `ray` and `mie`: `<ch>_time` (s since
  2000-01-01), `_lat`, `_lon`, `_z_bot`, `_z_top` (m), `_hlos`, `_hlos_err`
  (cm/s; missing: 32767 in `_hlos`, 65535 in `_hlos_err`), `_obs_type`,
  `_valid`, `_az`; global `Sources` (L2B files, baselines, processor
  versions).
- Profiles are rebuilt from the observations (the files have no profile
  identifier): see `obs/_limb_aeolus.py`. HLOS sign: `hlos_wind`'s
  `sign_convention` (found from ERA5, not ESA's documentation).
- The files cover 2018-09 to 2023-04, with no data from 2022-10-05 to
  2022-11-26.
- Checked: the 15th of every month 2020-02 to 2023-04, and 2018-12-31,
  2019-01-04 to 08, 2020-07-01, 2022-06-01; against ERA5 on three days. The separate
  `Aeolus/matlab/*_aeolus_hloswind.mat` files (2018-09 to 2022-01,
  Rayleigh only) are not read: they do not match these files observation by
  observation.
- Provenance: ESA Aeolus L2B (ALD_U_N_2B) through VirES (aeolus.services).

**Your own files** (`Misc`): any MAT-files in the ACE format above, named
`*<identifier>*yyyydddd*`, below `misc_path`.

### Nadir sounders (`load_nadir_sounders`)

**AIRS 3-D temperature retrievals** (Aqua)

- Files: daily `airs_yyyydddd.nc` (all 240 granules, dimension `time`),
  else granule files `airs_yyyy_ddd_ggg.nc`; usual layout `yyyy/` and
  `yyyy/ddd/`; netCDF.
- Variables: [along x cross] (135 x 90 per granule) `l1_lat`, `l1_lon`,
  `l1_time` (s since 2000-01-01T00:00Z); `ret_temp` (K) [along x cross x
  level]; `ret_z` (km, 27 levels). Fill -999 (padded granules), -9999
  positions, 0/0 positions, zero times; none declared.
- Checked: files of 2002-2025 on one machine; no version is stated in the
  files.
- Provenance: the retrieval of L. Hoffmann and M. J. Alexander (UNVERIFIED:
  reference to be added); placeholder, to be updated later: who made these
  files, from which data, and how to obtain them.

**AIRS 2-D brightness temperatures** (Aqua)

- Files: one per granule, `airs_yyyy_ddd_ggg.nc`; usual layout `yyyy/ddd/`
  or `ddd/`; netCDF.
- Variables: [NTRACK x NXTRACK] (135 x 90) `time` (s since
  2000-01-01T00:00Z), `lat`, `lon`, and for each channel group (`4mu`,
  `8mu`, `15mu_low`, `15mu_high`) `bt_<group>` (K), `bt_<group>_pt` (K) and
  `bt_<group>_var` (K^2) (no perturbation or variance for `8mu`).
- Checked: test copies only.
- Provenance: L. Hoffmann's AIRS brightness-temperature files; placeholder,
  to be updated later: where the originals are published and which version.

**IASI 3-D temperature retrievals** (MetOp)

- Files: one per orbit dump, `iasi_yyyymmddhhmmss.nc` (start of the first
  scan); usual layout `yyyy/ddd/`; only files whose name has the day's date
  are used for that day; netCDF.
- Variables: [rows x 60] `l1_lat`, `l1_lon`, `l1_time` (s since
  2000-01-01T00:00Z); `ret_temp` (K) [rows x 60 x level]; `ret_z` (km).
  Fills as for AIRS 3-D.
- Checked: files of 2007-2018 on one machine; no version is stated in the
  files.
- Provenance: the retrieval of L. Hoffmann (UNVERIFIED: reference to be
  added); placeholder, to be updated later: who made these files, from which
  data, and how to obtain them.

**SWANS model output sampled as AIRS**: AIRS 3-D granule files
(`SWANS/AIRS/dyamond_*/yyyy/ddd/airs_yyyy_ddd_ggg.nc`) load with
`file_path` and `sampled_from` (checked: GEOS 3 km, 2020-01-20).

### Trajectories (`load_trajectories`)

**DWD high-resolution radiosondes** (`DWDsondes`, folder `sondes`)

- Files: as unzipped from the DWD CDC (`.../radiosondes/high_resolution/
  historical/<yyyy>/`): `produkt_sec_aero_yyyymmdd_yyyymmdd_<station>.txt`
  (one per station and period, usually a year), with `Metadaten_Sekunde_Aero_<station>.txt`
  (launch and end times, sonde type) and `Metadaten_Geographie_<station>.txt`
  (station position and name) beside it.
- Variables: `;`-separated rows per sample: reference time, seconds since
  launch, position, geopotential height (gpm), pressure (hPa), temperature
  and dew point (degrees C), humidity (%), wind speed (m/s) and direction
  (degrees, from); -999 missing; exactly 0.00 positions are also missing.
- Checked: 1993-06-15, 2004-01-01, 2020-01-07, 2024-12-31 (up to 14
  stations); against ERA5 on five days (`docs/validation.md`; times are
  UTC); all 409 product files are in reference-time order. The field meanings are inferred from
  the files (no description is published with them): see
  `obs/_traj_dwd_sondes.py`.
- Provenance: DWD Climate Data Center, CC BY (as noted beside the data;
  credit DWD); `how_i_got_this_data.txt` beside the data.
- Two folders hold a copy of the same file: `1995/` has older, shorter
  copies of ten stations' 1993 files. The copy in the folder of its own
  year is used, else the largest (recorded in `multiple_matches`).

**IAGOS** (in-service aircraft, including MOZAIC; folder `IAGOS`)

- Files: one per flight, `IAGOS_timeseries_<yyyymmddHHMMSS><nn>[_L2_<v>].nc[4]`
  in `TimeSeries/iagos_<yyyy>/`, unzipped or in a zip archive (`iagos.zip`
  on eepc-0184: 67,321 flights 1994-2024, 13 GB, about 111 GB unzipped).
  Unzipped files are preferred, silently; the archive is the backup for
  flights not unzipped. Copies are matched by flight (departure time and
  number), so any unzipped copy replaces the archive's; of several
  unzipped copies the newest L2 version is used (in `multiple_matches`).
- Variables (every 4 s): `UTC_time` (s since the departure day), `lat`,
  `lon`, `gps_alt_AC`, `baro_alt_AC` (m), `air_press_AC` (Pa), `air_temp_AC`
  and the package's `air_temp_P1` (L2 3.1) or `air_temp_PM` (MOZAIC era)
  (K), `zon_wind_AC`, `mer_wind_AC` (m s-1), trace gases; validity flags
  `<var>_validity_flag` (0 good, 2 limited, 3 erroneous, 4 not validated, 7
  missing) in L2 3.1 files, `<var>_val` (all -9999 in the files checked) in
  older ones; fill -9999 (and -9999.9, e.g. in `baro_alt_AC`); one 2023
  flight has -999900 Pa (treated as missing). About 2% of L2 3.1 flights
  have a package temperature whose every value is flagged 3 or 7: their
  `air_temp_AC` is used. The archive holds 29 flights of January 2023 in
  both `iagos_2022/` and `iagos_2023/`; one copy is read.
- Checked: flights of 1998-02-27, 2010-12-04, 2019-03-18 and 2023-11-15
  (versions MOZAIC L2 revision 1.0, L2 3.1.0-3.1.2), against ERA5.
  MOZAIC-era files have no GPS altitude. One 2023-11-15 flight (D-AIKE) has
  82 take-off samples of exactly 173.16 K marked good; left as distributed.
- Provenance: IAGOS data portal (iagos-data.fr); the archive holds the
  version used by Wright et al. (2024, ACP; UNVERIFIED: reference not
  checked) and has not been updated since (`readme.txt` beside it).

**FAAM** (research aircraft; folder `FAAM`)

- Files: one per flight, `core_faam_<yyyymmdd>_v<nnn>_r<rev>_<flight>.nc`
  (FAAM core data, netCDF).
- Variables: `Time` (s since the flight day, 1 Hz); measurements (Time,
  sps<N>) at N per second with `<name>_FLAG` (bit masks, 0 good) and fill
  -9999: `LAT_GIN`, `LON_GIN`, `ALT_GIN` (m), `PS_RVSM` (hPa), `TAT_DI_R`,
  `TAT_ND_R` (K), `U_C`, `V_C`, `W_C` (m s-1), `PALT_RVS` (m) and about 150
  more (some at 2, 4, 20 or 64 per second). Of several revisions
  (`_r<rev>_`) of a flight, the highest is read.
- Checked: the 27 flights of 2025-07-09 to 2025-08-01 on eepc-0184 (v005,
  revision 0, processing software 25.0.0); three days against ERA5.
- Provenance: FAAM core data (the files state the UK Government Open
  Licence); placeholder, to be updated later: where these copies were
  downloaded from.

### Time-height series (`load_time_height`)

**ALOMAR lidar**: the same files as for `load_limb_sounders`, above.
