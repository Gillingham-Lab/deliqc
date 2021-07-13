import os


def get_rcparams():
    if os.path.exists("matplotlib.rcparams"):
        return "matplotlib.rcparams"
    else:
        return os.path.join(os.path.dirname(__file__), "..", "..", "res", "matplotlib.rcparams")