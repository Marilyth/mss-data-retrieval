# mss-data-retrieval
Scripts to get and process data for MSS on ecGates

# Scripts Setup
1. Clone this repository and move into it
```
git clone https://github.com/Marilyth/mss-data-retrieval.git
cd mss-data-retrieval
```
2. Make the shell scripts executable
```
chmod +x ./bin/*.sh
```
3. Install all requirements
```
pip --user -r install requirements.txt
```
4. Create a .bashrc for paths
```
cat .bashrc > ~/.bashrc
```

# Usage
1. After completing both setups, you can use this script as follows:
```
./bin/get_data.sh <date> <time>
```
For example, to get ERA5 data for March 2nd 2020 at 12 o'clock, use
```
./bin/get_data.sh 2020-03-02 12:00:00
```
