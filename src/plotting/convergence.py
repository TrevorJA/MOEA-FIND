"""Convergence diagnostics for MOEA-FIND optimization runs.

Tracks how the Pareto front evolves with function evaluations.
Uses a callback-based approach compatible with platypus and Borg.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt

from src.plotting.style import COLORS, apply_style


class ConvergenceTracker:
    """Tracks Pareto front evolution during optimization.

    Attach to an optimization run to record snapshots of the archive
    at regular intervals. Provides methods for plotting convergence.

    Args:
        snapshot_interval: Record a snapshot every N evaluations.
        objective_keys: Names of objectives (for labeling).
    """

    def __init__(
        self,
        snapshot_interval: int = 500,
        objective_keys: Optional[List[str]] = None,
    ):
        self.snapshot_interval = snapshot_interval
        self.objective_keys = objective_keys or []
        self.snapshots: List[Dict] = []
        self._eval_count = 0

    def record(self, nfe: int, archive_objs: np.ndarray):
        """Record a snapshot of the current archive.

        Args:
            nfe: Current number of function evaluations.
            archive_objs: Objective values of archive members, shape (n, k+1).
        """
        self.snapshots.append({
            "nfe": nfe,
            "n_solutions": len(archive_objs),
            "objectives": archive_objs.copy(),
        })

    @staticmethod
    def from_platypus_run(
        problem,
        epsilons: List[float],
        nfe_total: int,
        snapshot_interval: int = 500,
        seed: int = 42,
    ) -> Tuple["ConvergenceTracker", list]:
        """Run EpsNSGAII with convergence tracking.

        Runs the optimizer in chunks, recording snapshots at each interval.

        Args:
            problem: platypus Problem instance (with function set).
            epsilons: Epsilon values for EpsNSGAII.
            nfe_total: Total NFE budget.
            snapshot_interval: NFE between snapshots.
            seed: Random seed.

        Returns:
            Tuple of (tracker, final_result) where final_result is the
            list of Pareto solutions.
        """
        from platypus import EpsNSGAII
        import numpy as np

        np.random.seed(seed)
        tracker = ConvergenceTracker(snapshot_interval)
        algorithm = EpsNSGAII(problem, epsilons=epsilons)

        nfe_done = 0
        while nfe_done < nfe_total:
            chunk = min(snapshot_interval, nfe_total - nfe_done)
            algorithm.run(chunk)
            nfe_done += chunk

            if len(algorithm.result) > 0:
                objs = np.array([list(s.objectives) for s in algorithm.result])
                tracker.record(nfe_done, objs)

        return tracker, algorithm.result

    def plot_archive_size(
        self,
        ax: Optional[plt.Axes] = None,
        figsize: Tuple[float, float] = (8, 4),
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Plot number of Pareto solutions vs NFE.

        Args:
            ax: Optional axes to plot on.
            figsize: Figure size if creating new figure.

        Returns:
            (fig, ax) tuple.
        """
        apply_style()
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        nfes = [s["nfe"] for s in self.snapshots]
        sizes = [s["n_solutions"] for s in self.snapshots]

        ax.plot(nfes, sizes, "o-", color=COLORS["empirical"], markersize=3)
        ax.set_xlabel("Function Evaluations")
        ax.set_ylabel("Archive Size")
        ax.set_title("Pareto Archive Growth")
        return fig, ax

    def plot_objective_range_evolution(
        self,
        obj_indices: Optional[List[int]] = None,
        figsize: Tuple[float, float] = (12, 4),
    ) -> Tuple[plt.Figure, np.ndarray]:
        """Plot min/max range of each objective over NFE.

        Args:
            obj_indices: Which objectives to plot (default: all except Manhattan norm).
            figsize: Figure size.

        Returns:
            (fig, axes) tuple.
        """
        apply_style()
        if not self.snapshots:
            fig, ax = plt.subplots()
            return fig, np.array([ax])

        n_objs = self.snapshots[0]["objectives"].shape[1]
        if obj_indices is None:
            obj_indices = list(range(n_objs - 1))  # skip Manhattan norm

        fig, axes = plt.subplots(1, len(obj_indices), figsize=figsize)
        if len(obj_indices) == 1:
            axes = [axes]

        nfes = [s["nfe"] for s in self.snapshots]

        for ax, oi in zip(axes, obj_indices):
            mins = [s["objectives"][:, oi].min() for s in self.snapshots]
            maxs = [s["objectives"][:, oi].max() for s in self.snapshots]
            medians = [np.median(s["objectives"][:, oi]) for s in self.snapshots]

            ax.fill_between(nfes, mins, maxs, alpha=0.2, color=COLORS["empirical"])
            ax.plot(nfes, medians, "-", color=COLORS["empirical"], linewidth=1.5)
            ax.plot(nfes, mins, ":", color=COLORS["empirical"], linewidth=0.8)
            ax.plot(nfes, maxs, ":", color=COLORS["empirical"], linewidth=0.8)

            if self.objective_keys and oi < len(self.objective_keys):
                ax.set_title(self.objective_keys[oi])
            else:
                ax.set_title(f"Objective {oi+1}")
            ax.set_xlabel("NFE")
            ax.set_ylabel("Value")

        fig.suptitle("Objective Range Evolution", fontsize=12)
        fig.tight_layout()
        return fig, np.array(axes)

    def plot_front_snapshots(
        self,
        obj_x: int = 0,
        obj_y: int = 1,
        n_snapshots: int = 5,
        figsize: Tuple[float, float] = (8, 6),
    ) -> Tuple[plt.Figure, plt.Axes]:
        """Plot Pareto front at selected NFE snapshots.

        Args:
            obj_x: Objective index for x-axis.
            obj_y: Objective index for y-axis.
            n_snapshots: Number of snapshots to show.
            figsize: Figure size.

        Returns:
            (fig, ax) tuple.
        """
        apply_style()
        fig, ax = plt.subplots(figsize=figsize)

        if not self.snapshots:
            return fig, ax

        # Select evenly spaced snapshots
        total = len(self.snapshots)
        step = max(1, total // n_snapshots)
        selected = self.snapshots[::step]
        if self.snapshots[-1] not in selected:
            selected.append(self.snapshots[-1])

        cmap = plt.cm.viridis
        for i, snap in enumerate(selected):
            color = cmap(i / max(len(selected) - 1, 1))
            alpha = 0.3 + 0.7 * (i / max(len(selected) - 1, 1))
            objs = snap["objectives"]
            ax.scatter(objs[:, obj_x], objs[:, obj_y],
                       s=8, alpha=alpha, color=color,
                       label=f"NFE={snap['nfe']}")

        ax.set_xlabel(
            self.objective_keys[obj_x] if self.objective_keys and obj_x < len(self.objective_keys)
            else f"Objective {obj_x+1}"
        )
        ax.set_ylabel(
            self.objective_keys[obj_y] if self.objective_keys and obj_y < len(self.objective_keys)
            else f"Objective {obj_y+1}"
        )
        ax.set_title("Pareto Front Evolution")
        ax.legend(fontsize=8, loc="upper right")
        fig.tight_layout()
        return fig, ax
