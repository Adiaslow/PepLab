class CPKColors:
    """CPK coloring scheme for atoms."""

    COLORS = {
        'H': '#FFFFFF',   # White
        'He': '#D9FFFF',  # Light cyan
        'Li': '#CC80FF',  # Light purple
        'Be': '#C2FF00',  # Lime green
        'B': '#FFB5B5',   # Pink
        'C': '#909090',   # Grey
        'N': '#3050F8',   # Blue
        'O': '#FF0D0D',   # Red
        'F': '#90E050',   # Light green
        'Ne': '#B3E3F5',  # Cyan
        'Na': '#AB5CF2',  # Purple
        'Mg': '#8AFF00',  # Bright green
        'Al': '#BFA6A6',  # Light brown
        'Si': '#F0C8A0',  # Peach
        'P': '#FF8000',   # Orange
        'S': '#FFFF30',   # Yellow
        'Cl': '#1FF01F',  # Green
        'Ar': '#80D1E3',  # Light blue
        'K': '#8F40D4',   # Purple
        'Ca': '#3DFF00',  # Bright green
        'Sc': '#E6E6E6',  # Light grey
        'Ti': '#BFC2C7',  # Grey
        'V': '#A6A6AB',   # Grey
        'Cr': '#8A99C7',  # Grey-blue
        'Mn': '#9C7AC7',  # Purple-grey
        'Fe': '#E06633',  # Orange-brown
        'Co': '#F090A0',  # Pink
        'Ni': '#50D050',  # Green
        'Cu': '#C88033',  # Brown
        'Zn': '#7D80B0',  # Grey-blue
        'Ga': '#C28F8F',  # Pink-brown
        'Ge': '#668F8F',  # Grey-green
        'As': '#BD80E3',  # Purple
        'Se': '#FFA100',  # Orange
        'Br': '#A62929',  # Dark red
        'Kr': '#5CB8D1',  # Light blue
        'Rb': '#702EB0',  # Dark purple
        'Sr': '#00FF00',  # Green
        'Y': '#94FFFF',   # Light cyan
        'Zr': '#94E0E0',  # Light blue
        'Nb': '#73C2C9',  # Light blue
        'Mo': '#54B5B5',  # Blue-green
        'Tc': '#3B9E9E',  # Blue-green
        'Ru': '#248F8F',  # Blue-green
        'Rh': '#0A7D8C',  # Dark blue-green
        'Pd': '#006985',  # Dark blue
        'Ag': '#C0C0C0',  # Silver
        'Cd': '#FFD98F',  # Light yellow
        'In': '#A67573',  # Brown
        'Sn': '#668080',  # Grey-blue
        'Sb': '#9E63B5',  # Purple
        'Te': '#D47A00',  # Brown
        'I': '#940094',   # Purple
        'Xe': '#429EB0',  # Blue
        'Cs': '#57178F',  # Dark purple
        'Ba': '#00C900',  # Green
        'La': '#70D4FF',  # Light blue
        'Ce': '#FFFFC7',  # Very light yellow
        'Pr': '#D9FFC7',  # Very light green
        'Nd': '#C7FFC7',  # Light green
        'Pm': '#A3FFC7',  # Light green
        'Sm': '#8FFFC7',  # Light green
        'Eu': '#61FFC7',  # Light green
        'Gd': '#45FFC7',  # Light green
        'Tb': '#30FFC7',  # Light green
        'Dy': '#1FFFC7',  # Light green
        'Ho': '#00FF9C',  # Green
        'Er': '#00E675',  # Green
        'Tm': '#00D452',  # Green
        'Yb': '#00BF38',  # Green
        'Lu': '#00AB24',  # Green
        'Hf': '#4DC2FF',  # Light blue
        'Ta': '#4DA6FF',  # Light blue
        'W': '#2194D6',   # Blue
        'Re': '#267DAB',  # Blue
        'Os': '#266696',  # Blue
        'Ir': '#175487',  # Dark blue
        'Pt': '#D0D0E0',  # Grey
        'Au': '#FFD123',  # Gold
        'Hg': '#B8B8D0',  # Grey
        'Tl': '#A6544D',  # Brown
        'Pb': '#575961',  # Dark grey
        'Bi': '#9E4FB5',  # Purple
        'Po': '#AB5C00',  # Brown
        'At': '#754F45',  # Dark brown
        'Rn': '#428296',  # Blue
        'Fr': '#420066',  # Dark purple
        'Ra': '#007D00',  # Dark green
        'Ac': '#70ABFA',  # Light blue
        'Th': '#00BAFF',  # Light blue
        'Pa': '#00A1FF',  # Light blue
        'U': '#008FFF',   # Blue
        'Np': '#0080FF',  # Blue
        'Pu': '#006BFF',  # Blue
        'Am': '#545CF2',  # Blue
        'Cm': '#785CE3',  # Blue
        'Bk': '#8A4FE3',  # Purple
        'Cf': '#A136D4',  # Purple
        'Es': '#B31FD4',  # Purple
        'Fm': '#B31FBA',  # Purple
        'Md': '#B30DA6',  # Purple
        'No': '#BD0D87',  # Purple
        'Lr': '#C70066',  # Dark pink
        'Rf': '#CC0059',  # Dark pink
        'Db': '#D1004F',  # Dark pink
        'Sg': '#D90045',  # Dark pink
        'Bh': '#E00038',  # Dark pink
        'Hs': '#E6002E',  # Dark pink
        'Mt': '#EB0026',  # Dark pink
    }

    def get_color(self, cls, element: str) -> str:
        """Gets the CPK color for an element."""
        return cls.COLORS.get(element, '#808080')  # Default to grey
