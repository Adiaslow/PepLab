class CPKColors:
    """CPK coloring scheme for atoms."""

    COLORS = {
        'H': '#FFFFFF',  # White
        'C': '#909090',  # Grey
        'N': '#3050F8',  # Blue
        'O': '#FF0D0D',  # Red
        'F': '#90E050',  # Light green
        'Cl': '#1FF01F', # Green
        'Br': '#A62929', # Dark red
        'I': '#940094',  # Purple
        'P': '#FF8000',  # Orange
        'S': '#FFFF30',  # Yellow
        'B': '#FFB5B5',  # Pink
        'Na': '#0000FF', # Blue
        'K': '#8F40D4',  # Purple
        'Fe': '#FFA500', # Orange
    }

    @classmethod
    def get_color(cls, element: str) -> str:
        """Gets the CPK color for an element."""
        return cls.COLORS.get(element, '#808080')  # Default to grey
