# mss-data-retrieval
Scripts to get and process ERA5 data for MSS on Linux

# CDS-API Setup
1. Create an account at https://cds.climate.copernicus.eu/user/register \
Log into your account
2. Navigate to https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key \
Copy the content of the black box to the right and paste it into ~/.cdsapirc
5. Accept https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products

# Scripts Setup
1. Clone this repository and move into it
```
git clone https://github.com/Marilyth/mss-data-retrieval.git
cd mss-data-retrieval
```
2. Make the shell scripts executable
```
chmod +x ./bin/get_data.sh
chmod +x ./bin/add_pressure_gph.sh
```
3. Install all requirements
```
pip -r install requirements.txt
```

# Usage
After completing both setups, you can use this script as follows:
```
./bin/get_data.sh <date> <time>
```
For example, to get ERA5 data for March 2nd 2020 at 12 o'clock, use
```
./bin/get_data.sh 2020-03-02 12:00:00
```
Done, copy the .nc files to your mss data directory.
