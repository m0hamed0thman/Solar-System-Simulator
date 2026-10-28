# 🌌 Solar System Simulator

A parallel N-body simulation of the Solar System developed in C using **MPI (Message Passing Interface)**, with an interactive Python viewer for visualizing the motion of planets and other celestial bodies.

The simulator calculates gravitational interactions between celestial bodies, updates their positions and velocities over time, and saves the results to a CSV file for visualization.

## ✨ Features

- **N-body Simulation:** Simulates gravitational interactions between celestial bodies.
- **Parallel Computing:** Uses MPI to distribute computation across multiple processes.
- **Numerical Integration:** Uses the Velocity Verlet integration method to update positions and velocities.
- **Interactive Visualization:** Displays simulated trajectories using Python and Matplotlib.
- **CSV Data Processing:** Reads initial conditions from a CSV file and exports simulation results.
- **Customizable Input:** Supports different initial conditions through `initial.csv`.

## 🛠️ Technologies Used

- **C** — Simulation and gravitational calculations.
- **MPI / Open MPI** — Parallel processing.
- **Python 3** — Visualization.
- **NumPy** — Numerical data processing.
- **Matplotlib** — Interactive plotting and animation.
- **CSV** — Input and output data storage.

## 📋 Requirements

Before running the project, make sure you have installed:

- GCC
- Open MPI
- Python 3
- NumPy
- Matplotlib
- Git

### Install Dependencies on Ubuntu

Open a terminal and run:

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev python3 python3-pip python3-venv git
```

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/m0hamed0thman/Solar-System-Simulator.git
```

Navigate to the project directory:

```bash
cd Solar-System-Simulator
```

### 2. Set Up the Python Environment

Create a virtual environment:

```bash
python3 -m venv .venv
```

Activate it:

```bash
source .venv/bin/activate
```

Install the required Python libraries:

```bash
pip install numpy matplotlib
```

### 3. Compile the Simulation

Compile the MPI program using `mpicc`:

```bash
mpicc -O2 -Wall -Wextra -o planets_mpi planets_mpi_fixed.c -lm
```

This command generates the executable named `planets_mpi`.

**Note:** The repository also contains `planets_mpi.c`. The command above uses `planets_mpi_fixed.c`; compile only one implementation at a time.

### 4. Run the Simulation

Run the simulation using four MPI processes:

```bash
mpirun -np 4 ./planets_mpi initial.csv traj.csv
```

The program reads the initial conditions from `initial.csv` and writes the calculated trajectories to `traj.csv`.

You can change the number of processes by replacing `4` with another suitable number.

For example:

```bash
mpirun -np 2 ./planets_mpi initial.csv traj.csv
```

**Important:** The simulation runs for five simulated years by default. The computation may take some time, depending on your computer and the number of celestial bodies.

### 5. Visualize the Results

After the simulation finishes, run:

```bash
python3 sim_viewer.py traj.csv
```

An interactive visualization window will open, allowing you to explore the simulated Solar System.

If you want to display calendar dates instead of simulation days, specify the starting date:

```bash
python3 sim_viewer.py traj.csv --epoch 2026-01-01
```

The epoch specifies the calendar date corresponding to simulation time `t = 0`.

## 🎮 Visualization Controls

Use the following controls to interact with the simulation:

| Control | Action |
|---|---|
| Mouse hover | Display information about a celestial body |
| Left-click | Focus the view on a celestial body |
| Mouse wheel | Zoom in or out |
| `Space` | Play or pause the animation |
| `Left Arrow` | Move to the previous frame |
| `Right Arrow` | Move to the next frame |
| `+` / `-` | Zoom in or out |
| `Home` | Return to the full-system view |
| `T` | Toggle trajectory trails |
| `L` | Toggle body labels |
| `S` | Save a snapshot |
| `F` | Toggle fullscreen |
| `Q` | Quit the viewer |

## 📁 Project Structure

```text
Solar-System-Simulator/
│
├── initial.csv             # Initial conditions of celestial bodies
├── traj.csv                # Generated simulation trajectories
├── fetch_all_solar_system.py # Script for retrieving Solar System data
├── planets_mpi.c            # MPI-based simulation implementation
├── planets_mpi_fixed.c      # Updated simulation implementation
├── planets_mpi              # Compiled executable, if present
├── sim_viewer.py             # Interactive Python visualization
├── nbody_planets.html        # HTML-based visualization or related resource
├── README.md                 # Project documentation
└── .vscode/                  # VS Code configuration
```

Some files, including `traj.csv` and the compiled executable, may be generated locally rather than tracked in Git.

## 📊 Input and Output

### Input: `initial.csv`

The input file contains the initial physical parameters of the simulated celestial bodies.

The simulation supports these formats:

**Nine-column format:**

```csv
name,parent,gm_km3_s2,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud
```

**Ten-column format:**

```csv
name,parent,gm_km3_s2,gm_au3_d2,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud
```

Where:

- `name` — Name of the celestial body.
- `parent` — Parent body, where applicable.
- `gm_km3_s2` — Gravitational parameter in km³/s².
- `gm_au3_d2` — Gravitational parameter in AU³/day²; accepted in the ten-column format.
- `x_au`, `y_au`, `z_au` — Initial position coordinates in astronomical units (AU).
- `vx_aud`, `vy_aud`, `vz_aud` — Initial velocity components in AU/day.

The simulation converts the gravitational parameter into AU³/day² for its internal calculations.

**Note:** The column names and order must match the supported format.

### Output: `traj.csv`

The simulation generates a CSV file containing the calculated position and velocity of each body at recorded simulation times.

The output format is:

```csv
t_days,name,x_au,y_au,z_au,vx_aud,vy_aud,vz_aud
```

Where:

- `t_days` — Simulation time in days.
- `name` — Name of the celestial body.
- `x_au`, `y_au`, `z_au` — Position coordinates in AU.
- `vx_aud`, `vy_aud`, `vz_aud` — Velocity components in AU/day.

The generated file is used by `sim_viewer.py` to animate the simulated motion.

## ⚙️ How It Works

The project follows these main steps:

1. **Load Initial Conditions:** Read the celestial bodies and their initial positions, velocities, and gravitational parameters from `initial.csv`.
2. **Initialize MPI:** Start multiple processes to distribute the simulation workload.
3. **Calculate Gravitational Acceleration:** Calculate the gravitational acceleration acting on each assigned body due to the other bodies.
4. **Update Positions and Velocities:** Apply the Velocity Verlet integration method to advance the system through time.
5. **Synchronize Processes:** Use MPI communication to share updated simulation data between processes.
6. **Save Results:** Write the recorded positions and velocities to `traj.csv`.
7. **Visualize the Simulation:** Load the trajectory data into the Python viewer to explore the motion interactively.

## 🧪 Troubleshooting

### `mpicc: command not found`

Install Open MPI:

```bash
sudo apt install openmpi-bin libopenmpi-dev
```

### `No module named 'numpy'` or `No module named 'matplotlib'`

Activate your virtual environment and install the dependencies:

```bash
source .venv/bin/activate
pip install numpy matplotlib
```

### `cannot open initial.csv`

Make sure you are running the simulation from the project directory and that `initial.csv` exists.

You can also provide the input file's path explicitly:

```bash
mpirun -np 4 ./planets_mpi /path/to/initial.csv traj.csv
```

### The visualization cannot find `traj.csv`

Run the simulation first, then launch the viewer:

```bash
mpirun -np 4 ./planets_mpi initial.csv traj.csv
python3 sim_viewer.py traj.csv
```

### Open MPI reports that there are not enough slots

Try using fewer processes:

```bash
mpirun -np 2 ./planets_mpi initial.csv traj.csv
```

Alternatively, for a local educational experiment, you can allow Open MPI to oversubscribe the available slots:

```bash
mpirun --oversubscribe -np 4 ./planets_mpi initial.csv traj.csv
```

Oversubscription allows more MPI processes than the available slots; it does not guarantee better performance.

## ⚠️ Performance Notes

This project is primarily an educational demonstration of parallel computing and gravitational simulation.

Because the current dataset contains a relatively small number of celestial bodies, MPI communication overhead may outweigh the benefits of parallelization. Running the program with more processes does not necessarily make it faster.

## 🎯 Project Goals

- Explore N-body gravitational simulation.
- Practice parallel programming with MPI.
- Understand numerical integration using Velocity Verlet.
- Process scientific data using CSV files.
- Visualize celestial motion using Python and Matplotlib.

## 📌 Future Improvements

Potential improvements include:

- Add more celestial bodies and customizable simulation parameters.
- Compare serial and parallel execution times.
- Measure speedup and parallel efficiency for different process counts.
- Improve numerical accuracy and stability.
- Add more visualization options and simulation controls.
- Allow users to configure the simulation duration and time step.

## 👨‍💻 Author

**Mohamed Hassan Othman**

GitHub: [@m0hamed0thman](https://github.com/m0hamed0thman)

## 📄 License

No license is currently specified in this repository. Check the repository for an existing license before redistributing or reusing the project.