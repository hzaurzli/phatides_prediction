import os
import argparse


def suffix_rename(path_input,suffix):
		path = path_input
		if path[-1] == '/':
			path = path
		else:
			path = path + '/'

		files = os.listdir(path)
		for filename in files:
			portion = os.path.splitext(filename)
			if portion[1] == suffix:
				newname = portion[0] + ".fna"
				os.rename(path + filename,path + newname)


if __name__ == '__main__':
		parser = argparse.ArgumentParser(description="suffix rename")
		parser.add_argument("-p", "--path", required=True, type=str, help="genome sequence path")
		parser.add_argument("-s", "--suffix", required=True, type=str, help="old suffix/To be modified suffix")
		Args = parser.parse_args()

		path = Args.path
		suffix = Args.suffix

		if suffix[0] == '.':
			suffix = suffix
		else:
			suffix = '.' + suffix

		suffix_rename(path,suffix)
