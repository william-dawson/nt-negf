from sys import argv
from matplotlib import pyplot as plt

if __name__ == "__main__":
    x = []
    y = []
    with open(argv[1]) as ifile:
        for line in ifile:
            split = [float(x) for x in line.split()]
            x.append(27.2114 * split[0])
            y.append(split[1])

    fig, axs = plt.subplots(1, 1, figsize=(6, 3))
    axs.plot(x, y, 'o--')
    axs.set_ylim(-.1, 1.1 * max(y))
    axs.set_xlabel("Energy (eV)", fontsize=12)
    axs.set_ylabel("Transmission", fontsize=12)
    plt.tight_layout()
    plt.show()
