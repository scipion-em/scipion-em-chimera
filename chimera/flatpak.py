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

def ask_install_scope():
    while True:
        answer = input(
            "Installation type:\n"
            "  1) User (--user)\n"
            "  2) System-wide (all users, needs sudo)\n"
            "Choose [1/2]: "
        ).strip()

        if answer == "1":
            return "user"
        elif answer == "2":
            return "system"

        print("Invalid choice. Please enter 1 or 2.\n")
