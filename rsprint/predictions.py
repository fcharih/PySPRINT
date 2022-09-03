#!/usr/bin/env python3

"""Module containing the Predictions class and related methods

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""
import pickle
from typing import Dict, List, Tuple, Union

import numpy as np

class Predictions(object):
    def __init__(self, scores: np.array, protein_names: List[str], peptide_names: List[str] = []) -> "Predictions":
        self.names = [x for x in protein_names] + [x for x in peptide_names]
        self.name_index_map = { p: i for i, p in enumerate(protein_names + peptide_names) }
        self.index_name_map = { i: p for i, p in enumerate(protein_names + peptide_names) }
        self.scores = scores

    @staticmethod
    def from_file(filename: str) -> "Predictions":
        with open(f"{filename}", "rb") as pickled_file:
            data = pickle.load(pickled_file)
            predictions = Predictions(data["scores"], data["names"])
            return predictions

    def get_scores(self, protein: str, format: str = "dict") -> Union[np.array, Dict]:
        scores_array = self.scores[self.name_index_map[protein]]
        if format == "array":
            return scores_array
        elif format == "dict":
            scores_dict = dict(zip(self.names, list(scores_array)))
            return scores_dict
        else:
            raise Exception(f"`{format}` is not a valid format. Please choose one from `dict` and `array`.")

    def get_score(self, protein1: str, protein2: str) -> float:
        protein1_index = self.name_index_map[protein1]
        protein2_index = self.name_index_map[protein2]
        score = self.scores[protein1_index, protein2_index]
        return score

    def save(self, filename: str):
        with open(filename, "wb") as pickle_file:
            file_content = {
                "names": self.names,
                "scores": self.scores
            }
            pickle.dump(file_content, pickle_file)