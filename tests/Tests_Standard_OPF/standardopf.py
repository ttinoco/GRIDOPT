import pfnet as pf
import gridopt as gopt
import csv
from gridopt.power_flow.method_error import PFmethodError
import timeit
listobjf=[]
listiterations=[]
listiterationsc=[]
listPowermis=[]
listTime=[]
listcase=['case9','case14','case39','case89pegase','case118','case300','case1354pegase','case2869pegase','case3375wp']
listmethod =['IpoptOPF','AugLOPF']
for case in listcase:
	i=0
	with open('sol_wo_thlim'+case+'.csv') as f:
	    reader = csv.reader(f)
	    for row in reader:
	    	
	    	if i==0:
	    		objfmatlab		=row[0]
	    	if i==1:
	    		iterationmatlab = row[0]
	    	if i==2:
	    		maxPmismatlab 	= row[0]
	    	if i==3:
	    		timematlab 		= row[0]
	    	i +=1
	net = pf.Network()
	pfcase = 'pf'+case+'.mat'
	net.load(pfcase)
	listTimec=[]
	listobjfc=[]
	listiterationsc=[]
	listPowermisc=[]
	for methodname in listmethod:

		method = gopt.power_flow.new_method(methodname)
		try:	

			start = timeit.default_timer()

			method.set_parameters({'quiet': True})

			method.solve(net)

			stop = timeit.default_timer()

			results = method.get_results()


			Iterations     = results['iterations']
			Objfunction    = results['net properties']['gen_P_cost']
			Power_mismatch = results['net properties']['bus_P_mis']
			time 		   = stop-start

			listobjfc.append((100*(Objfunction-float(objfmatlab)))/float(objfmatlab))
			listiterationsc.append(Iterations)
			listPowermisc.append(Power_mismatch/net.base_power)
			listTimec.append(time)
		except PFmethodError:

			Iterations     = 0
			Objfunction    = 0
			Power_mismatch = 0
			time 		   = 0
			listobjfc.append((100*(Objfunction-float(objfmatlab)))/float(objfmatlab))
			listiterationsc.append(Iterations)
			listPowermisc.append(Power_mismatch/net.base_power)
			listTimec.append(time)


	listobjfc.append(objfmatlab)
	listiterationsc.append(iterationmatlab)
	listPowermisc.append(maxPmismatlab)	
	listTimec.append(timematlab)

	listobjf.append(listobjfc)
	listiterations.append(listiterationsc)
	listPowermis.append(listPowermisc)
	listTime.append(listTimec)

print('type of case','	----	IPOPT/AUGLOPF/MATLAB Objfunction','	----	IPOT/AUGLOPF/MATLAB iterations','	----	IPOPT/AUGLOPF/MATALB Time','	----	IPOPT/AUGLOPF/MATALB maxPowermis')
for i in range(0,len(listiterations)):
	print(listcase[i],'	----	',listobjf[i],'	----	',listiterations[i],'	----	',listTime[i],'	----	',listPowermis[i])
