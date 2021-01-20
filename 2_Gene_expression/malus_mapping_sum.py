ctl={}
for i in range(0,15):
	print x[i]
	y=open(x[i].strip()+'ReadsPerGene.out.tab').readlines()
	y=[j.split()[0] for j in y[5:] if int(j.split()[1])>100]
	ctl[x[i]]=y
	
ctl_sum=[]
for i in range(0,15):
	ctl_sum=list(set().union(ctl_sum,ctl[x[i]]))

#get significantly expressed gene in the experimental set and non-overlapping with the control
exp={}
for i in range(15,len(x)):
	print x[i]
	try:
		y=open(x[i].strip()+'ReadsPerGene.out.tab').readlines()
		y=[j.split()[0] for j in y[5:] if int(j.split()[1])>100 and not j.split()[0] in ctl_sum]
		exp[x[i]]=y
	except IOError:pass

exp_genes={}
for i in range(15,len(x)):
	try:
		for j in exp[x[i]]:
			try:
				exp_genes[j]=exp_genes[j]+1
			except KeyError:
				exp_genes[j]=1
	except KeyError:pass

for k in exp_genes.keys():
	if exp_genes[k]>2:
		print(k+'\t'+`exp_genes[k]`)