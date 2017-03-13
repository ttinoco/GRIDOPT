import pfnet as pf
import gridopt as gopt
import csv
listobjf=[]
listiterations=[]
listobjfc=[]
listiterationsc=[]
listPowermis=[]
listcase=['case9','case14','case39','case89pegase','case118','case300']#,'case1354pegase','case2869pegase','case3375wp']
for case in listcase:
	i=0
	with open('sol'+case+'.csv') as f:
	    reader = csv.reader(f)
	    for row in reader:
	    	
	    	if i==0:
	    		objfmatlab=row[0]
	    	if i==1:
	    		iterationmatlab = row[0]
	    	if i==2:
	    		maxPmismatlab = row[0]
	    	i +=1
	net = pf.Network()
	pfcase = 'pf'+case+'.mat'
	net.load(pfcase)

	method = gopt.power_flow.new_method('AugLOPF')
	method.solve(net)
	results = method.get_results()

	iterations_aug = results['iterations']
	objfunction_aug = results['net properties']['gen_P_cost']
	Power_mismatch_aug = results['net properties']['bus_P_mis']



	method1 = gopt.power_flow.new_method('IpoptOPF')
	method1.solve(net)
	results1 = method1.get_results()

	iterations_ipopt = results1['iterations']
	objfunction_ipopt = results1['net properties']['gen_P_cost']
	Power_mismatch_ipopt = results1['net properties']['bus_P_mis']
	


	listobjfc=[objfunction_ipopt,objfunction_aug, objfmatlab]
	listiterationsc=[iterations_ipopt,iterations_aug,iterationmatlab]
	listobjf.append(listobjfc)
	listiterations.append(listiterationsc)
	listPowermis.append([Power_mismatch_ipopt,Power_mismatch_aug, maxPmismatlab])
print('type of case','IPOPT/AUGLOPF/MATLAB Objfunction','IPOT/AUGLOPF/MATLAB iterations','IPOPT/AUGLOPF/MATALB maxPowermis')
for i in range(0,len(listiterations)):
	print(listcase[i],listobjf[i],listiterations[i],listPowermis[i])
"""
method.update_network(net)
base = net.base_power
print(net.get_load(2).P*base+net.get_load(1).P*base+net.get_load(0).P*base)


method = gopt.power_flow.new_method('IpoptOPF')
method.solve(net)
results = method.get_results()

iteration = results['iterations']
objfunction = results['net properties']['gen_P_cost']

print(objfunction, iteration)
"""