#!/usr/bin/env python3
"""
Plot SCF convergence data from scf_data.txt
"""

import sys
import matplotlib.pyplot as plt
import numpy as np


def plot_scf_data(filename="scf_data.txt"):
    """Plot mu vs trace convergence."""
    # Read data
    data = np.loadtxt(filename, skiprows=1)

    if data.size == 0:
        print(f"ERROR: No data found in {filename}")
        sys.exit(1)

    # Extract columns
    if data.ndim == 1:  # Single iteration
        iterations = np.array([data[0]])
        mu_history = np.array([data[1]])
        trace_history = np.array([data[2]])
    else:
        iterations = data[:, 0]
        mu_history = data[:, 1]
        trace_history = data[:, 2]

    # Compute target trace (should be constant, use first value as reference)
    # For display, we'll compute it from the data trend
    target_trace = np.mean([trace_history[0], trace_history[-1]])

    # Find the closest point to target
    errors = np.abs(trace_history - target_trace)
    final_idx = np.argmin(errors)
    final_mu = mu_history[final_idx]
    final_trace = trace_history[final_idx]
    final_error = errors[final_idx]

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # Plot 1: mu vs trace
    ax1.plot(
        mu_history,
        trace_history,
        "o-",
        label="SCF iterations",
        markersize=8,
        linewidth=2,
        alpha=0.7,
        color='tab:blue'
    )

    # Mark start and end
    ax1.plot(
        mu_history[0], trace_history[0], "go", markersize=12, label="Start", zorder=5
    )
    ax1.plot(
        mu_history[-1], trace_history[-1], "ro", markersize=12, label="Final", zorder=5
    )

    # Plot target trace
    ax1.axhline(
        y=target_trace,
        color="r",
        linestyle="--",
        linewidth=2,
        alpha=0.5,
        label=f"Target trace ≈ {target_trace:.4f}",
    )

    # Mark final mu
    ax1.axvline(
        x=final_mu,
        color="g",
        linestyle="--",
        linewidth=2,
        alpha=0.5,
        label=f"Final μ = {final_mu:.6f} Ha",
    )

    ax1.set_xlabel("Chemical Potential μ (Ha)", fontsize=12)
    ax1.set_ylabel("Trace", fontsize=12)
    ax1.set_title(
        "Self-Consistent Chemical Potential Search", fontsize=14, fontweight="bold"
    )
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)

    # Plot 2: Convergence history (iteration vs error)
    ax2.semilogy(
        iterations,
        errors,
        "o-",
        markersize=8,
        linewidth=2,
        alpha=0.7,
        color='tab:orange'
    )

    ax2.set_xlabel("Iteration", fontsize=12)
    ax2.set_ylabel("Trace Error (log scale)", fontsize=12)
    ax2.set_title("Convergence History", fontsize=14, fontweight="bold")
    ax2.grid(True, alpha=0.3, which='both')

    # Add final error annotation
    ax2.annotate(
        f"Final error: {final_error:.6f}",
        xy=(iterations[-1], errors[-1]),
        xytext=(10, 10),
        textcoords='offset points',
        fontsize=10,
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
    )

    plt.tight_layout()

    # Print summary
    print("\n" + "="*50)
    print("SCF Convergence Summary")
    print("="*50)
    print(f"Total iterations: {len(iterations)}")
    print(f"Target trace:     {target_trace:.6f}")
    print(f"Final μ:          {final_mu:.8f} Ha")
    print(f"Final trace:      {final_trace:.6f}")
    print(f"Final error:      {final_error:.6f}")
    print(f"μ range:          [{mu_history.min():.6f}, {mu_history.max():.6f}] Ha")
    print("="*50 + "\n")

    return fig


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Plot SCF convergence data")
    parser.add_argument(
        "filename",
        nargs="?",
        default="scf_data.txt",
        help="SCF data file (default: scf_data.txt)",
    )
    parser.add_argument(
        "--save",
        type=str,
        default=None,
        help="Save plot to file instead of displaying (e.g., scf_plot.png)",
    )
    args = parser.parse_args()

    try:
        fig = plot_scf_data(args.filename)

        if args.save:
            fig.savefig(args.save, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {args.save}")
        else:
            print("Displaying plot... (close window to exit)")
            plt.show()

    except FileNotFoundError:
        print(f"ERROR: File '{args.filename}' not found")
        print("Make sure you've run the SCF driver first to generate scf_data.txt")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
