import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
	file = sys.argv[1]
else:
	file = 'results/benchmark.csv'


"""
print specific variants of the data
"""

data = pd.read_csv(file, sep=';')
print(data["Algorithm"].unique())
print(data["Variant"].unique())

algorithms = ["Multiplication2"]
variants = ["Multiplication2",'Multiplication2_0','Multiplication2_1',
 'Multiplication2_2']
colors = ["green", "blue", "orange", "red", "gray"]
line_styles = ["-", "-", "-", "-", "-"]
labels = [ "Baseline", "Opt: 0", "Opt: 1", "Opt: 2", "truncatedMul (CAMPARY)"]
markers = ["o", "v","^", "s","D"]
baseline_performance = data.query('Algorithm == "Multiplication2"').query('Variant == "Multiplication2"')["Performance"].iloc[-1]

input_sizes = []
i = 0
for algo in algorithms:
	algo_query = 'Algorithm == "%s"' % algo
	algo_data = data.query(algo_query)
	for variant in variants:
		print("Algorithm: {}\nVariant: {}".format(algo, variant))
		variant_data = data.query('Variant == "%s"' % variant)
		if len(variant_data) < 1:
			raise Exception("Could not find variant: {}".format(variant))
		plt.plot(variant_data['Input size'], 
	   				variant_data['Performance'], 
					label=variant, 
					color = colors[i], 
					linestyle = line_styles[i], 
	   				marker = markers[i])
		if len(algo_data['Input size']) > len(input_sizes):
			input_sizes = algo_data['Input size']
		

		labels[i] = labels[i] + " (+{:.0f}%)".format((variant_data['Performance'].iloc[-1]/baseline_performance-1)*100) 	
		i+=1
plt.legend(loc='lower right', labels = labels)

plt.xticks(input_sizes, input_sizes)

#plt.title('Intel(R) Core(TM) i7-1065G7 CPU @ 1.3GHz\nFlags: -mfma -O3 -fno-tree-vectorize', size = 7)
#plt.suptitle('Performance of mult2', size = 15)

plt.ylabel('Performance [flops/cycle]')
plt.xlabel('Input Size')
plt.grid(axis = 'y')
plt.savefig('mult2_optimization.pdf')
plt.ylim(0)
plt.show()
