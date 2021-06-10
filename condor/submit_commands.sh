#!/bin/bash
-z 3 -L 50 -w 50  -N 10000 -b 25 -o HGCALTile_Gamma1p33MeV.root
-d 1.6 -r 3 -W 1.4 -R 0.4 -S 0.3 -a 1.0 -b 25 -w 50 -N 100

# Geometry stuff

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 30  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 20  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 50  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 30  -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 20 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d $(seq 1.5 0.1 2.0) -r $(seq 3.0 0.5 6.0) -L 50  -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.5 0.5 6.0)  -L 30 -W 3.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.5 0.1 2.0) -r $(seq 3.5 0.5 6.0)  -L 20 -W 3.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d $(seq 1.5 0.1 2.0) -r $(seq 3.5 0.5 6.0)  -L 50 -W 3.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985


# Reflectivity stuff
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.5 -r 3.0 -L 30  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d 1.5 -r 3.0 -L 20  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d 1.5 -r 3.0 -L 50  -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.5 -r 3.0 -L 30  -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d 1.5 -r 3.0 -L 20  -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d 1.5 -r 3.0 -L 50  -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 0.95 0.92 0.9 -p 0.2 0.5 0.7 0.9

# Radiation damage stuff
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.5 -r 3.0 -L 30  -W 1.4 -R 0.4 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d 1.5 -r 3.0 -L 20  -W 1.4 -R 0.4 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d 1.5 -r 3.0 -L 50  -W 1.4 -R 0.4 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.5 -r 3.0 -L 30  -W 2.0 -R 0.1 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d 1.5 -r 3.0 -L 20  -W 2.0 -R 0.1 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0  24 ) -d 1.5 -r 3.0 -L 50  -W 2.0 -R 0.1 -S 0.3 -a 1.0 0.9 0.8 0.5 0.1 -w 0.985

# New optimization commands
python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 2.0 -R 0.1 -S 0.3 -a 0.20 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 0.20 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 2.0 -R 0.1 -S 0.3 -a 0.15 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 0.15 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 2.0 -R 0.1 -S 0.3 -a 0.10 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 0.10 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 1.4 -R 0.4 -S 0.3 -a 1.0 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 1.4 -R 0.4 -S 0.3 -a 0.20 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 1.4 -R 0.4 -S 0.3 -a 0.20 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 1.4 -R 0.4 -S 0.3 -a 0.15 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 1.4 -R 0.4 -S 0.3 -a 0.15 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0  9.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 20 -W 1.4 -R 0.4 -S 0.3 -a 0.10 -w 0.985
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d $(seq 1.6 0.2 2.0) -r $(seq 6.0 0.2 7.0)  -L 30 -W 1.4 -R 0.4 -S 0.3 -a 0.10 -w 0.985

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.8 -r 6.0  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 -A $(seq 0.01 0.01 0.05)
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.8 -r 6.0  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 -D $(seq 0.10 0.10 0.50)

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.8 -r 6.0  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 0.15 -w 0.985 -A $(seq 0.01 0.01 0.05)
python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 1.8 -r 6.0  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 0.15 -w 0.985 -D $(seq 0.10 0.10 0.50)

python3 SubmitJobs.py -x $(seq 0 1.0 14.5) -d 0.9 -r 6.0  -L 30 -W 2.0 -R 0.1 -S 0.3 -a 1.0 -w 0.985 -T 0 1 2 3
