"""
hplc_sim.optimization
======================
Optional extension: search for HPLC conditions that maximize overall
resolution / minimize run time, using either Monte Carlo random search
or a lightweight genetic algorithm (no external GA library dependency).
"""
from __future__ import annotations
import random
import copy
from typing import Dict, List, Tuple

from .simulator import HPLCSimulator


def _score(result: Dict, target_rs: float = 1.5) -> float:
    """Fitness: reward resolved peaks, penalize long run time and under-resolved pairs."""
    peaks = result["peaks"]
    if not peaks:
        return -1e9
    rs_values = [p.get("Rs_vs_prev") for p in peaks if p.get("Rs_vs_prev") is not None]
    if not rs_values:
        rs_penalty = 0.0
    else:
        rs_penalty = sum(max(0.0, target_rs - rs) ** 2 for rs in rs_values)
    run_time = result["method"]["run_time_min"]
    return -rs_penalty - 0.02 * run_time  # higher is better


def monte_carlo_optimize(
    sim: HPLCSimulator,
    n_iter: int = 200,
    phi_bounds: Tuple[float, float] = (0.05, 0.95),
    flow_bounds: Tuple[float, float] = (0.3, 2.0),
) -> Dict:
    """Random search over isocratic phi and flow rate."""
    best = None
    best_score = -1e18
    for _ in range(n_iter):
        phi = random.uniform(*phi_bounds)
        flow = random.uniform(*flow_bounds)
        sim.patch_method(isocratic_phi=phi, mode="isocratic")
        sim.patch_flow(flow)
        result = sim.run(add_noise=False, n_points=1500)
        s = _score(result)
        if s > best_score:
            best_score = s
            best = {"phi": phi, "flow_ml_min": flow, "score": s, "result": result}
    return best


def genetic_optimize(
    sim: HPLCSimulator,
    population_size: int = 24,
    generations: int = 15,
    phi_bounds: Tuple[float, float] = (0.05, 0.95),
    flow_bounds: Tuple[float, float] = (0.3, 2.0),
    mutation_rate: float = 0.2,
) -> Dict:
    """Simple real-valued GA (tournament selection, blend crossover, gaussian mutation)."""
    def random_individual():
        return [random.uniform(*phi_bounds), random.uniform(*flow_bounds)]

    def evaluate(ind):
        sim.patch_method(isocratic_phi=ind[0], mode="isocratic")
        sim.patch_flow(ind[1])
        result = sim.run(add_noise=False, n_points=1200)
        return _score(result), result

    population = [random_individual() for _ in range(population_size)]
    best_overall = None
    best_overall_score = -1e18

    for gen in range(generations):
        scored = [(ind, *evaluate(ind)) for ind in population]
        scored.sort(key=lambda x: x[1], reverse=True)
        if scored[0][1] > best_overall_score:
            best_overall_score = scored[0][1]
            best_overall = {"phi": scored[0][0][0], "flow_ml_min": scored[0][0][1],
                             "score": scored[0][1], "result": scored[0][2]}

        # elitism: keep top 20%
        n_elite = max(1, population_size // 5)
        new_population = [ind for ind, _, _ in scored[:n_elite]]

        # tournament selection + blend crossover
        while len(new_population) < population_size:
            p1 = min(random.sample(scored, 3), key=lambda x: -x[1])[0]
            p2 = min(random.sample(scored, 3), key=lambda x: -x[1])[0]
            alpha = random.random()
            child = [alpha * p1[i] + (1 - alpha) * p2[i] for i in range(2)]
            if random.random() < mutation_rate:
                child[0] += random.gauss(0, 0.05)
                child[1] += random.gauss(0, 0.1)
            child[0] = min(max(child[0], phi_bounds[0]), phi_bounds[1])
            child[1] = min(max(child[1], flow_bounds[0]), flow_bounds[1])
            new_population.append(child)

        population = new_population

    return best_overall
