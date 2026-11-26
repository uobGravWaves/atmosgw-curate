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

- the daily-file readers (HIRDLS, MLS, SABER, SOFIE, MIPAS) search the whole
  folder for a name containing the date, as their MATLAB originals do;
- the others (ACE, GNSS, planetary-wave files, AIRS, IASI) try the usual path
  first (e.g. `yyyy/ace_v52_yyyydddd.mat`, `yyyy/ddd/airs_yyyy_ddd_ggg.nc`),
  then search for that name anywhere below the folder (one folder listing per
  load: under a second for the largest folders checked, AIRS and IASI, about
  50 000 entries each);
- ALOMAR files are looked for in the folder itself; Sampler files in
  `sample_path` itself.

If more than one file has the name or date looked for, the first in search
order is used (at each level, the contents of sub-folders before the files
beside them, each in name order) and the
choice is warned and listed in the output attribute `multiple_matches`; the
newest MLS version is chosen instead of the first. A file at its usual path is
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
- Checked: files of 2007-2016 on one machine; no version is stated in the
  files.
- Provenance: the retrieval of L. Hoffmann (UNVERIFIED: reference to be
  added); placeholder, to be updated later: who made these files, from which
  data, and how to obtain them.

### Time-height series (`load_time_height`)

**ALOMAR lidar**: the same files as for `load_limb_sounders`, above.
