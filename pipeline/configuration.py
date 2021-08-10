from collections import OrderedDict


def get_indelible_config():
	indelible_config = OrderedDict()

	indelible_config["[TYPE]"] = 'AMINOACID 2'
	indelible_config["[MODEL]"] = 'modelname'
	indelible_config["[submodel]"] = 'WAG'
	indelible_config["[indelmodel]"] = 'POW 1.7 500'
	indelible_config["[indelrate]"] = '0.0'
	indelible_config["[rates]"] = ' 0.25 0.50 10'
	indelible_config["[statefreq]"] = ' 0.25  0.25  0.25  0.25'
	indelible_config["[TREE]"] = 'treename (A:0.1,B:0.1);'
	indelible_config["[PARTITIONS]"] = 'partitionname\n[treename modelname 3000]'
	indelible_config["[EVOLVE]"] = "partitionname 100 outputname" + " " # note: indlible requires space at last command.

	return indelible_config

def get_sparta_config():
	sparta_config = OrderedDict()

	# sparta_config["_indelibleTemplateControlFile"] = "control_indelible_template.txt"
	# sparta_config["_dawgTemplateControlFile"] = ""
	# sparta_config["_dawgSimulator"] = "0"
	sparta_config["_inputRealMSAFile"] = "results/ref_msa.aa.fasta"
	sparta_config["_inputTreeFileName"] = "results/RAxML_tree.tree"
	sparta_config["_outputGoodParamsFile"] = "results/SpartaABC_data_name_model_name.posterior_params"
	sparta_config["_numberOfSamplesToKeep"] = "100"
	sparta_config["_alignmentMode"] = "0"
	sparta_config["_similarity_mode"] = "0"
	sparta_config["_modelType"] = "eq"
	sparta_config["_minRLVal"] = "50"
	sparta_config["_maxRLVal"] = "500.0"
	sparta_config["_minIRVal"] = "0"
	sparta_config["_maxIRVal"] = "0.05"
	sparta_config["_distanceCutOff"] = "1.0"

	sparta_config["_wAvgUniqueGapSize"] = "-1.0"
	sparta_config["_wMSAMin"] = "-1.0"
	sparta_config["_wNumGapsLenTwo"] = "-1.0"
	sparta_config["_wAvgGapSize"] = "-1.0"
	sparta_config["_wTotNumGaps"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFour"] = "-1.0"
	sparta_config["_wNumGapsLenOne"] = "-1.0"
	sparta_config["_wMSAMax"] = "-1.0"
	sparta_config["_wMSALen"] = "-1.0"
	sparta_config["_wTotNumUniqueGaps"] = "-1.0"
	sparta_config["_wNumGapsLenThree"] = "-1.0"
	sparta_config["_wNumGapsLenOneIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenOneIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenOneInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenTwoInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenThreeInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourIn1Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourIn2Pos"] = "-1.0"
	sparta_config["_wNumGapsLenAtLeastFourInNMinus1Pos"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_0_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_1_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_2_gaps"] = "-1.0"
	sparta_config["_wNumberOfMSA_position_with_n_minus_1_gaps"] = "-1.0"

	sparta_config["_only_real_stats"] = "0"
	sparta_config["_alignments_output"] = "10"
	sparta_config["_outputAlignmnetsFile"] = "results/alignments.fasta"

	sparta_config["_numSimulations"] = "100000"
	sparta_config["_numBurnIn"] = "1"



	return sparta_config

# configuration classes should be used as singleton, so do not change anything inside after instanciation in main.
# TODO: change to hidden members and add getters.
class general_config:

	def __init__(self,pipeline_path, results_path, tree_file_name, msa_file_name, available_models, skip_config, clean_run, verbose, op_sys):
		self.pipeline_path = pipeline_path
		self.results_path = results_path
		self.msa_file_name = msa_file_name
		self.tree_file_name = tree_file_name
		self.available_models = available_models
		self.skip_config = skip_config
		self.clean_run = clean_run
		self.verbose = verbose
		self.op_sys = op_sys

class simulations_config:
	def __init__(self, simulations_num, burnin_num ,min_ir,max_ir, alignments_num):
		self.simulations_num = simulations_num
		self.burnin_num = burnin_num
		self.min_ir = min_ir
		self.max_ir = max_ir
		self.alignments_num = alignments_num


class correction_config:
	def __init__(self, model_parameters, filter_p, regression_mode):
		self.model_parameters = model_parameters
		self.filter_p = filter_p
		self.regression_mode = regression_mode


class inference_config:
	def __init__(self, lib, number_top, size_threshold):
		self.lib = lib
		self.number_top = number_top
		self.size_threshold = size_threshold
		
		self.advanced = False # Enable neural network and regression when True.


if __name__ == "__main__":
	sparta_config = get_sparta_config()
	sparta_config["_inputRealMSAFile"] = "/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/results/msa.fasta"
	sparta_config["_inputTreeFileName"] = "/home/elyawy/development/Msc/Thesis/Working_dir/SpartaABC/results/RAxML_tree.tree"
	with open("test.conf",'w') as f:
		for key in get_sparta_config():
			to_write = f'{key} {sparta_config[key]}\n'
			f.write(to_write)

