import dalmatian
import re
import pandas as pd

class FC:
	def __init__(self):
		self.ws_list = dalmatian.firecloud.api.list_workspaces().json();

	def _get_TCGA_workspace_names(self):
		ws_names = [(y["name"],
		             re.match(r"TCGA_([A-Z]+)_ControlledAccess.*", y["name"]), 
		             y["namespace"])
					for y in [x["workspace"] for x in self.ws_list]];
		ws_names = [(x[0], x[1].groups()[0], x[2]) for x in ws_names if x[1] != None];

		return pd.DataFrame(ws_names, columns = ["workspace", "ttype", "namespace"])

	def _get_all_workspace_names(self):
		return pd.DataFrame(
		  [(y["name"], y["namespace"]) for y in
		    [x["workspace"] for x in self.ws_list ]],
	      columns = ["workspace", "namespace"]
		)

fc = FC();
get_TCGA_workspace_names = fc._get_TCGA_workspace_names
get_all_workspace_names = fc._get_all_workspace_names
