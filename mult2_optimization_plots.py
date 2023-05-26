import pandas as pd
import matplotlib.pyplot as plt


data = pd.read_csv('results/benchmark.csv', sep=';')

algorithms = data['Algorithm'].unique()

algo = algorithms[1]
input_sizes = []
for algo in algorithms:
	algo_query = 'Algorithm == "%s"' % algo
	algo_data = data.query(algo_query)
	plt.plot(algo_data['Input size'], algo_data['Performance'], label=algo)
	if len(algo_data['Input size']) > len(input_sizes):
		input_sizes = algo_data['Input size']
plt.legend()
plt.xscale('log')
plt.xticks(input_sizes, input_sizes)
plt.title('Performance [flops/cycle] of mult2')
plt.xlabel('Input Size')
plt.savefig('mult2_optimization.png')
plt.show()
