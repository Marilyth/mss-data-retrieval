# mss-data-retrieval
Scripts to get and process data for MSS on ecGates

# Scripts Setup
1. Clone this repository and move into it
```
git clone https://github.com/Marilyth/mss-data-retrieval.git
cd mss-data-retrieval
git checkout ecgate-forecast
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
source ~/.bashrc
```

# Usage
1. After completing both setups, you can use this script as follows:
```
./bin/get_data.sh <date> <time> <step>
```
Where \<date\> and \<time\> is the date and time where the forecast was created, and \<step\> is how many hours after the date and time you want your data.\
For example, to get a forecast created at 22nd of April 2021 at 0 o'clock, for 23rd and 24th of April at 0 o'clock, use
```
./bin/get_data.sh 2021-04-22 00:00:00 24/48
```
