# PepLab

PepLab is a peptide library generation toolkit designed to support various workflows for processing peptide libraries, generating 3D structures, analyzing reactivity, planning reactions, and executing reactions. It provides a comprehensive platform for peptide-related research and analysis.

## Features

- **Design**: Create and enumerate peptide libraries using a variety of design strategies including:
  - Combinatorial methods
  - Generative AI models
  - Genetic algorithms
  - Markov chain Monte Carlo simulations
- **Analysis**: Analyze library data inlcuding:
  - Peptide properties
  - Experimental results
  - DEL sequencing data
  - ...
- **Modeling**: Model peptide structures and interactions using:
  - State-of-the-art structure prediction models
  - Molecular dynamics
  - Molecular docking
  - ...
- **Optimization**: Optimize peptide libraries for desired properties using:
  - Genetic algorithms
  - Machine learning models
  - Markov chain Monte Carlo simulations

## Installation

1. Clone the repository:

```bash
git clone https://github.com/yourusername/peplab.git
cd peplab
```

2. Create and activate a virtual environment:

On Unix/macOS

```bash
python -m venv venv
source venv/bin/activate
```

On Windows

```bash
python -m venv venv
.\venv\Scripts\activate
```

3. Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

1. Activate the virtual environment:

On Unix/macOS:

```bash
source venv/bin/activate
```

On Windows (untested):

```bash
.\venv\Scripts\activate
```

2. Run the application:

```bash
flask run
```

The application will be available at http://127.0.0.1:5000

To stop the server:

1. Press CTRL+C in the terminal
2. Deactivate the virtual environment:

```bash
deactivate
```

## Development

To run in development mode with auto-reload:

```bash
flask run --debug
```

### Contributing

Contributions are welcome! Here are the steps:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make your changes.
4. Commit your changes (`git commit -am 'Add new feature'`).
5. Push to the branch (`git push origin feature-branch`).
6. Create a new Pull Request.

### License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

### Acknowledgements

- Andre
