"""
Natural Selection Simulation
----------------------------
Simulation of natural selection in a giraffe population.
Each giraffe has an 8-base DNA sequence:
  - bases 0-3 -> neck score (can it reach the food?)
  - bases 4-7 -> lifespan score

Usage:
    python simulation.py            # runs the 3 experiments, saves figures in figures/
"""

import os
import random

import matplotlib
matplotlib.use("Agg")  # no window needed, figures are saved to files
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# 1. DNA and constants
# ---------------------------------------------------------------------------
A, T, G, C = "A", "T", "G", "C"
bases = [A, T, G, C]
complement = {A: T, T: A, G: C, C: G}

m_prob = 0.08  # mutation probability per base (changed by the experiments)

valeur_cou = {A: 3, T: 1, G: 4, C: 0}   # value of each base for the neck score
valeur_life = {A: 3, T: 3, G: 2, C: 4}  # value of each base for the lifespan

AGE_MIN = 4            # minimum age for reproduction
INTERVALLE_REPRO = 3   # years between two births for a female
sexe = ["male", "female"]


# ---------------------------------------------------------------------------
# 2. DNA functions
# ---------------------------------------------------------------------------
def replication_mut(adn, m_prob):
    """Copy a DNA sequence; each base can mutate with probability m_prob."""
    adn_replique = []
    for base in adn:
        if random.random() < m_prob:
            adn_replique.append(random.choice([b for b in bases if b != base]))
        else:
            adn_replique.append(base)
    return adn_replique


def score_eat(adn):
    """Neck score = sum of the values of the first 4 bases."""
    allele = adn[0:4]
    return sum(valeur_cou[base] for base in allele)


def score_life(adn):
    """Lifespan score = sum of the values of the last 4 bases."""
    allele = adn[4:8]
    return sum(valeur_life[base] for base in allele)


def reproduction(g1, g2):
    """Return the DNA of a child (None if both parents have the same sex)."""
    if g1.sexe == g2.sexe:
        return None
    a1 = replication_mut(g1.adn, m_prob)
    a2 = replication_mut(g2.adn, m_prob)
    gene_cou = random.choice([a1[0:4], a2[0:4]])
    gene_vie = random.choice([a1[4:8], a2[4:8]])
    return gene_cou + gene_vie


# ---------------------------------------------------------------------------
# 3. Classes
# ---------------------------------------------------------------------------
class Girafe:
    def __init__(self, adn, sexe):
        self.adn = adn
        self.sexe = sexe
        self.score_cou = score_eat(adn)
        self.score_vie = score_life(adn)
        self.age = 0
        self.vivant = True
        self.cooldown = 0

    def vieillir(self, affame=False):
        """Age by 1 year (2 if hungry). Dies when age >= lifespan score."""
        self.age += 2 if affame else 1
        if self.cooldown > 0:
            self.cooldown -= 1
        if self.age >= self.score_vie:
            self.vivant = False


class Environnement:
    def __init__(self, high):
        self.high_food = high  # height of the food (tree height)

    def girafe_mange(self, girafe):
        return girafe.score_cou >= self.high_food


def trouver_partenaire(girafe, population):
    """Pick a random living adult of the opposite sex (None if none)."""
    candidats = [g for g in population
                 if g.vivant and g.sexe != girafe.sexe and g.age >= AGE_MIN]
    return random.choice(candidats) if candidats else None


def creer_girafe_aleatoire():
    """Giraffe with a random DNA sequence and a random sex."""
    adn = [random.choice(bases) for _ in range(8)]
    return Girafe(adn, random.choice(sexe))


def creer_clone(i):
    """Identical giraffe (neck score 5, long lifespan), alternating sex."""
    adn = [A, T, T, C, C, C, C, C]
    return Girafe(adn.copy(), "male" if i % 2 == 0 else "female")


# ---------------------------------------------------------------------------
# 4. Simulation
# ---------------------------------------------------------------------------
def simulation_stat(population, tour, environnement, prob_repro, verbose=False):
    """
    Run the simulation for `tour` years.
    Returns the final population and a history (one dict per year).
    """
    historique = []
    for a in range(tour):
        newborn = []
        for girafe in population:
            if not girafe.vivant:
                continue
            girafe.vieillir(affame=not environnement.girafe_mange(girafe))
            if not girafe.vivant:
                continue
            if (girafe.sexe == "female" and girafe.age >= AGE_MIN
                    and girafe.cooldown == 0 and random.random() < prob_repro):
                partenaire = trouver_partenaire(girafe, population)
                if partenaire is not None:
                    enfant_adn = reproduction(girafe, partenaire)
                    if enfant_adn is not None:
                        enfant_sexe = random.choice(sexe)
                        newborn.append(Girafe(enfant_adn, enfant_sexe))
                        girafe.cooldown = INTERVALLE_REPRO
        population = population + newborn
        population = [g for g in population if g.vivant]
        historique.append({
            "tour": a,
            "n": len(population),
            "scores_cou": [g.score_cou for g in population],
            "genes": ["".join(g.adn[0:4]) for g in population],
            "scores_life": [g.score_vie for g in population],
        })
        if verbose:
            print(f"Year {a}: {len(population)} giraffes alive")
    return population, historique


def lancer(high, m, repetitions=5, tours=50):
    """
    Run `repetitions` random populations with tree height `high` and mutation
    rate `m`. Returns the mean neck score and the mean number of gene variants,
    year by year (averaged over the repetitions).
    """
    global m_prob
    m_prob = m
    scores, divs = [], []
    for _ in range(repetitions):
        pop = [creer_girafe_aleatoire() for _ in range(50)]
        env = Environnement(high=high)
        pop, hist = simulation_stat(pop, tours, env, 0.65)
        scores.append([np.mean(h["scores_cou"]) if h["n"] > 0 else np.nan
                       for h in hist])
        divs.append([len(set(h["genes"])) if h["n"] > 0 else np.nan
                     for h in hist])
    return np.nanmean(scores, axis=0), np.nanmean(divs, axis=0)


# ---------------------------------------------------------------------------
# 5. Experiments (one per section of the paper)
# ---------------------------------------------------------------------------
def experience_cou(outdir):
    """3.1 - Neck size evolution: 15 random populations, high = 9."""
    global m_prob
    m_prob = 0.08
    moyenne_graph = {}
    for i in range(15):
        env = Environnement(high=9)
        pop = [creer_girafe_aleatoire() for _ in range(50)]
        pop, hist = simulation_stat(pop, 50, env, 0.65)
        moyenne_graph[i] = [np.mean(h["scores_cou"]) if h["n"] > 0 else np.nan
                            for h in hist]

    plt.figure(figsize=(10, 6))
    cmap = plt.get_cmap("tab20")
    for i, courbe in moyenne_graph.items():
        plt.plot(courbe, color=cmap(i), label=f"Simulation {i + 1}")
    plt.xlabel("Year")
    plt.ylabel("Mean neck score")
    plt.title("Mean neck score evolution (15 simulations, high = 9)")
    plt.legend(ncol=3, fontsize=8)
    plt.savefig(os.path.join(outdir, "fig1_neck_score.png"), dpi=300)
    plt.close()


def experience_mutation(outdir):
    """3.2 - Mutation rate: identical low-neck giraffes, high = 9."""
    global m_prob
    mutations = [0, 0.001, 0.01, 0.05, 0.1, 0.3]
    repetitions = 5
    res_cou, res_pop = {}, {}

    for m in mutations:
        m_prob = m
        cou_reps, pop_reps = [], []
        for _ in range(repetitions):
            pop = [creer_clone(i) for i in range(50)]
            env = Environnement(high=9)
            pop, hist = simulation_stat(pop, 50, env, 0.65)
            cou_reps.append([np.mean(h["scores_cou"]) if h["n"] > 0 else np.nan
                             for h in hist])
            pop_reps.append([h["n"] for h in hist])
        res_cou[m] = np.nanmean(cou_reps, axis=0)
        res_pop[m] = np.mean(pop_reps, axis=0)
        survivantes = sum(p[-1] > 0 for p in pop_reps)
        print(f"m = {m}: {survivantes}/{repetitions} populations survived")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for m in mutations:
        axes[0].plot(res_cou[m], label=f"m = {m}")
        axes[1].plot(res_pop[m], label=f"m = {m}")
    axes[0].set_title("Mean neck score")
    axes[0].set_ylabel("Neck score")
    axes[1].set_title("Population size")
    axes[1].set_ylabel("Number of giraffes")
    for ax in axes:
        ax.set_xlabel("Year")
    axes[0].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig2_mutation_rate.png"), dpi=300)
    plt.close()


def experience_interaction(outdir):
    """3.3 - Interaction between environmental pressure and mutation rate."""
    mutations = [0, 0.001, 0.01, 0.05, 0.1, 0.3]
    resultats = {}
    for high in [0, 9]:
        for m in mutations:
            resultats[(high, m)] = lancer(high, m)
            print(f"high = {high}, m = {m}: done")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for col, high in enumerate([0, 9]):
        for m in mutations:
            sc, dv = resultats[(high, m)]
            axes[0, col].plot(sc, label=f"m = {m}")
            axes[1, col].plot(dv, label=f"m = {m}")
        axes[0, col].set_title(f"high = {high}")
        axes[0, col].set_ylabel("Mean neck score")
        axes[1, col].set_ylabel("Number of variants of gene 0-4")
        axes[1, col].set_xlabel("Year")
    axes[0, 0].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "fig3_interaction.png"), dpi=300)
    plt.close()


# ---------------------------------------------------------------------------
# 6. Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    random.seed(42)  # remove this line for different results each run
    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)

    print("Experiment 1: neck size...")
    experience_cou(outdir)
    print("Experiment 2: mutation rate...")
    experience_mutation(outdir)
    print("Experiment 3: interaction...")
    experience_interaction(outdir)
    print(f"Done. Figures saved in ./{outdir}/")
