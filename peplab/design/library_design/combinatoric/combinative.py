<<<<<<< HEAD
from itertools import product
from typing import List, Callable, Any, Optional, Dict

class Combinative:
    def __init__(self, *item_lists: List[List[Any]]):
        """
        Initializes the Combinative class with lists of items to be combined.

        Parameters:
        - item_lists: A variable number of lists, each containing items to combine.
        """
        self.item_lists = item_lists

    def generate_combinations(self, filter_func: Optional[Callable[[List[Any]], bool]] = None) -> List[List[Any]]:
        """
        Generates all possible combinations of items from the input lists.

        Parameters:
        - filter_func: Optional function that takes a combination (a list) as input and returns True if the 
                       combination should be included in the output, False otherwise.

        Returns:
        - List of lists, where each sublist is a valid combination of items.
        """
        # Generate all combinations using Cartesian product
        all_combinations = list(product(*self.item_lists))
        
        # Apply filter if provided
        if filter_func:
            filtered_combinations = [list(comb) for comb in all_combinations if filter_func(comb)]
            return filtered_combinations
        
        # Return all combinations without filtering
        return [list(comb) for comb in all_combinations]

    def to_dict(self, combinations: List[List[Any]]) -> Dict[int, List[Any]]:
        """
        Converts a list of combinations into a dictionary where keys are indices and values are combinations.

        Parameters:
        - combinations: List of combinations to convert.

        Returns:
        - Dictionary with index keys and combination values.
        """
        return {i: comb for i, comb in enumerate(combinations)}

    '''def save_to_file(self, combinations: List[List[Any]], filename: str) -> None:
        """
        Saves combinations to a text file, each combination on a new line.

        Parameters:
        - combinations: List of combinations to save.
        - filename: Name of the file to save combinations to.
        """
        with open(filename, "w") as file:
            for comb in combinations:
                line = ", ".join(map(str, comb))
                file.write(f"{line}\n")
        print(f"Combinations saved to {filename}")'''

# Example Usage
if __name__ == "__main__":
    # Initialize with three lists of items
    comb = Combinative([1, 2], ['a', 'b'], ['X', 'Y'])

    # Generate all combinations
    all_combs = comb.generate_combinations()
    print("All Combinations:", all_combs)

    # Generate combinations with a filter (e.g., only combinations where first item is 1)
    filtered_combs = comb.generate_combinations(filter_func=lambda x: x[0] == 1)
    print("Filtered Combinations:", filtered_combs)

    # Convert to dictionary
    comb_dict = comb.to_dict(all_combs)
    print("Combinations Dictionary:", comb_dict)

    # Save combinations to a file
    #comb.save_to_file(all_combs, "combinations.txt")
=======

>>>>>>> origin/main
