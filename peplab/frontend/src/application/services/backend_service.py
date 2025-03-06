import requests

BASE_URL = "http://localhost:5000/api/combinatorial"

def fetch_permutations(sequence, num_permutations):
    response = requests.post(f"{BASE_URL}/permutations", json={
        "sequence": sequence,
        "num_permutations": num_permutations
    })
    return response.json()