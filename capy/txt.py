import pandas as pd
import os
import typing
import subprocess

def parsein(X, col, regex, fields):
	T = parse(X[col], regex, fields)
	return pd.concat([X, T], 1)

def parse(X, regex, fields):
	T = X.str.extract(regex).rename(columns = dict(enumerate(fields)))
	return T

def print_full(D):
	with pd.option_context('display.max_rows', None, 'display.max_columns', None,
	  'display.max_colwidth', 1000):
		print(D)

def to_clipboard(string: str):
	p = subprocess.Popen("xclip", stdin = subprocess.PIPE)
	p.communicate(str.encode(string))
