""" Auxiliary functions to handle flatpak installation and check if the application is installed. """
import subprocess


def is_installed(app_id):
    result = subprocess.run(
        ["flatpak", "info", "--user", app_id],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return result.returncode == 0


def install_flatpak(package):
    subprocess.run(
        ["flatpak", "install", "--user", "-y", "--noninteractive", package],
        check=True
    )