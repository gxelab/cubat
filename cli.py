import click


@click.group()
@click.version_option()
def cli():
    """Naval Fate.
    This is the docopt example adopted to Click but with some actual
    commands implemented and not just the empty parsing which really
    is not all that interesting.
    """


@cli.group()
def ship():
    """Manages ships."""
