import argh
from deliqc.cli import about, run, helpers, plot

commands = [
    about.about,
    run.run,
    plot.plot,
]


def main():
    p = argh.ArghParser()
    p.add_commands(sorted(commands, key=lambda x: x.__name__))

    p.set_default_command(about.about)

    try:
        argh.dispatch(p)
    except helpers.Error as e:
        print(e)
        return 1
