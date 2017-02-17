### Author: Edward Huang

import os

### This script generates the beginning directories.

def main():
    for folder in ['./data', './results']:
        if not os.path.exists(folder):
            os.makedirs(folder)

if __name__ == '__main__':
    main()