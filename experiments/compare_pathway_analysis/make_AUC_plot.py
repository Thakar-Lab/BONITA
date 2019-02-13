import csv
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log
import seaborn as sns

# complete and save a plot
def finishPlot(xlabeler, ylabeler, plotname):
	plt.ylabel(ylabeler, labelpad=0)
	plt.xlabel(xlabeler, labelpad=1.)
	sns.despine()
	sns.set_style("ticks")
	sns.set(font_scale=1) 
	fonter=8
	plt.rc('font', size=fonter)          # controls default text sizes
	plt.rc('axes', titlesize=fonter)     # fontsize of the axes title
	plt.rc('axes', labelsize=fonter)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=fonter)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=fonter)    # fontsize of the tick labels
	plt.rc('legend', fontsize=fonter)    # legend fontsize
	plt.rc('figure', titlesize=fonter)  # fontsize of the figure title
	plt.savefig(plotname, bbox_inches='tight', transparent=True, dpi=600)
	plt.savefig(plotname[:-3]+'svg', bbox_inches='tight', transparent=True, dpi=600)
	plt.clf()

def styleSetter():
	sns.set_context("paper")
	sns.set_style("ticks")
	sns.set(font_scale=1) 
	fonter=8
	plt.rc('font', size=fonter)          # controls default text sizes
	plt.rc('axes', titlesize=fonter)     # fontsize of the axes title
	plt.rc('axes', labelsize=fonter)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=fonter)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=fonter)    # fontsize of the tick labels
	plt.rc('legend', fontsize=fonter)    # legend fontsize
	plt.rc('figure', titlesize=fonter)  # fontsize of the figure title
	sns.set_style("ticks")

def genROC(tcomb,fcomb):
	labels=[1 for i in tcomb]
	labels2=[0 for i in fcomb]
	tcomber=list(tcomb)
	tcomber.extend(fcomb)
	labels.extend(labels2)
	
	fpr, tpr, scores = roc_curve(labels, tcomber)
	lw = 2
	return(fpr,tpr)


def rocPlotAll(bon, clip, cam, name):
	b_fpr,b_tpr=genROC(bon[1],bon[0])
	clip_fpr,clip_tpr=genROC(clip[1],clip[0])
	cam_fpr, cam_tpr=genROC(cam[1],cam[0])
	roc_auc_b = auc(b_fpr,b_tpr)
	roc_auc_clip = auc(clip_fpr,clip_tpr)
	roc_auc_cam = auc(cam_fpr, cam_tpr)
	lw=2
	styleSetter()
	colors=	sns.color_palette("colorblind", 3)
	plt.plot(b_fpr, b_tpr, color=colors[0], lw=lw, label='BONITA (area = %0.2f)' % roc_auc_b)
	plt.plot(clip_fpr, clip_tpr, color=colors[1], lw=lw, label='CLIPPER (area = %0.2f)' % roc_auc_clip)
	plt.plot(cam_fpr, cam_tpr, color=colors[2], lw=lw, label='CAMERA (area = %0.2f)' % roc_auc_cam)
	plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.legend(loc="lower right")
	finishPlot('False Positive Rate', 'True Positive Rate', name)

def rocPlotAll2(bon, clip, cam,bon2, clip2, cam2, name, label1, label2):
	b_fpr,b_tpr=genROC(bon[1],bon[0])
	clip_fpr,clip_tpr=genROC(clip[1],clip[0])
	cam_fpr, cam_tpr=genROC(cam[1],cam[0])
	roc_auc_b = auc(b_fpr,b_tpr)
	roc_auc_clip = auc(clip_fpr,clip_tpr)
	roc_auc_cam = auc(cam_fpr, cam_tpr)
	
	b_fpr2,b_tpr2=genROC(bon2[1],bon2[0])
	clip_fpr2,clip_tpr2=genROC(clip2[1],clip2[0])
	cam_fpr2, cam_tpr2=genROC(cam2[1],cam2[0])
	roc_auc_b2 = auc(b_fpr2,b_tpr2)
	roc_auc_clip2 = auc(clip_fpr2,clip_tpr2)
	roc_auc_cam2 = auc(cam_fpr2, cam_tpr2)
	lw=1.5
	styleSetter()
	plt.figure(figsize=(2.6,2))
	colors=	sns.color_palette("colorblind", 3)

	dot=mpl.lines.Line2D([0,0], [0,0], linewidth=lw, linestyle='dotted', color='black')
	solid=mpl.lines.Line2D([0,0], [0,0], linewidth=lw, linestyle='solid', color='black')
	
	plt.plot(cam_fpr, cam_tpr, color=colors[0], lw=lw, linestyle='dotted', label=label1+'CAMERA (%0.2f)' % roc_auc_cam)
	plt.plot(clip_fpr, clip_tpr, color=colors[1], lw=lw, linestyle='dotted', label=label1+'CLIPPER (%0.2f)' % roc_auc_clip)
	plt.plot(b_fpr, b_tpr, color=colors[2], lw=lw, linestyle='dotted', label=label1+'BONITA (%0.2f)' % roc_auc_b)

	plt.plot(cam_fpr2, cam_tpr2, color=colors[0], lw=lw,label=label2+'CAMERA (%0.2f)' % roc_auc_cam2)
	plt.plot(clip_fpr2, clip_tpr2, color=colors[1], lw=lw, label=label2+'CLIPPER (%0.2f)' % roc_auc_clip2)
	plt.plot(b_fpr2, b_tpr2, color=colors[2], lw=lw, label=label2+'BONITA (%0.2f)' % roc_auc_b2)
	print(roc_auc_b)
	print(roc_auc_cam)
	print(roc_auc_clip)
	plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.legend([dot,solid],['0.5','2.0'], framealpha=0.0,title="Induced (log2) \nAttenuation")
	finishPlot('False Positive Rate', 'True Positive Rate', name)



def ROCplotter(tcomb,fcomb,name):
	fpr,tpr=genROC(tcomb,fcomb)
	lw = 2
	roc_auc = auc(fpr, tpr)
	print(roc_auc)
	plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.savefig(name)
	plt.clf()



def readBONITA(filename):
	results=[]
	for i in [0,5,10,15,20]:
		results.append([])
	with open(filename, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in spamreader:
			for i in range(1,len(row)):
				if len(row[i])>0:
					results[i-1].append(float(row[i]))
	for q in results:
		q.pop(0)
	return(results)
def analyzeBONITA():
	# BONITA
	results=readBONITA("BONITA/NP.csv")
	# ROCplotter(results[1],results[0],'BONITA_NP_5.png')
	# ROCplotter(results[2],results[0],'BONITA_NP_10.png')
	# ROCplotter(results[3],results[0],'BONITA_NP_15.png')
	results2=readBONITA("BONITA/RS.csv")
	# ROCplotter(results2[1],results2[0],'BONITA_RS_5.png')
	# ROCplotter(results2[2],results2[0],'BONITA_RS_10.png')
	# ROCplotter(results2[3],results2[0],'BONITA_RS_15.png')
	return([results,results2])
def readClipperFile(filename):
	results=[]
	with open(filename, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
		for row in spamreader:
			results.append(row)
		results.pop(0)
		results.pop(0)
	results2=[]
	for lister in results:
		lister.pop(0)
		temper=[]
		for element in lister:
			if float(element)==0:
				temper.append(1000.)
			else:
				temper.append(-1.*log(float(element),10))
		results2.append(temper)
	return(results2)


def analyzeClipper():
	graphs=['hsa04060','hsa04612','hsa04660','hsa04650','hsa04630','hsa04350']
	topDir="clipper/"
	gen='NP'
	finalResult=[]
	for gen in ['NP','RS']:
		for result in ['comb','mean']:
			data=[[] for a in [0,5,10,15,20]]
			for graph in graphs:
				filename=topDir+graph+'_'+gen+'_'+result+'.csv'
				output=readClipperFile(filename)
				for i in range(len(output)):
					data[i].extend(output[i])
			# print(gen)
			# print(result)
			if result=='mean':
				finalResult.append(list(data))
			# ROCplotter(data[1],data[0],'CLIPPER_'+gen+'_'+result+'_5.png')
			# ROCplotter(data[2],data[0],'CLIPPER_'+gen+'_'+result+'_10.png')
			# ROCplotter(data[3],data[0],'CLIPPER_'+gen+'_'+result+'_15.png')
			# ROCplotter(data[4],data[0],'CLIPPER_'+gen+'_'+result+'_20.png')
	return(finalResult)

def analyzeCamera():
	graphs=['hsa04060','hsa04612','hsa04660','hsa04650','hsa04630','hsa04350']
	topDir="CAMERA/"
	gen='NP'
	datas=[]
	for gen in ['NP','RS']:
		data=[[] for a in [0,5,10,15,20]]
		for graph in graphs:
			filename=topDir+graph+'_'+gen+'_'+'CAMERA.csv'
			output=readClipperFile(filename)
			for i in range(len(output)):
				data[i].extend(output[i])
		# print(gen)
		datas.append(list(data))
		# ROCplotter(data[1],data[0],'CAMERA_'+gen+'_5.png')
		# ROCplotter(data[2],data[0],'CAMERA_'+gen+'_10.png')
		# ROCplotter(data[3],data[0],'CAMERA_'+gen+'_15.png')
		# ROCplotter(data[4],data[0],'CAMERA_'+gen+'_20.png')
	return(datas)

if __name__ == '__main__':
	graphs=['hsa04630']
	bon=analyzeBONITA()
	clip=analyzeClipper()
	cam=analyzeCamera()
	# rocPlotAll2([bon[0][0],bon[0][1]], [clip[0][0],clip[0][1]], [cam[0][0],cam[0][1]], [bon[0][0],bon[0][3]], [clip[0][0],clip[0][3]], [cam[0][0],cam[0][3]],'NP_all_5_15.png', '0.5 ', '1.5 ')
	# rocPlotAll2([bon[0][0],bon[0][1]], [clip[0][0],clip[0][1]], [cam[0][0],cam[0][1]], [bon[0][0],bon[0][2]], [clip[0][0],clip[0][2]], [cam[0][0],cam[0][2]],'NP_all_5_10.png', '0.5 ', '1.0 ')
	rocPlotAll2([bon[0][0],bon[0][1]], [clip[0][0],clip[0][1]], [cam[0][0],cam[0][1]], [bon[0][0],bon[0][4]], [clip[0][0],clip[0][4]], [cam[0][0],cam[0][4]],'NP_all_5_20_plos.png', '0.5 ', '2.0 ')
	# rocPlotAll2([bon[1][0],bon[1][1]], [clip[1][0],clip[1][1]], [cam[1][0],cam[1][1]], [bon[1][0],bon[1][3]], [clip[1][0],clip[1][3]], [cam[1][0],cam[1][3]],'RS_all_5_15.png', '0.5 ', '1.5 ')
	# rocPlotAll2([bon[1][0],bon[1][1]], [clip[1][0],clip[1][1]], [cam[1][0],cam[1][1]], [bon[1][1],bon[1][2]], [clip[1][0],clip[1][2]], [cam[1][0],cam[1][2]],'RS_all_5_10.png', '0.5 ', '1.0 ')
	rocPlotAll2([bon[1][0],bon[1][1]], [clip[1][0],clip[1][1]], [cam[1][0],cam[1][1]], [bon[1][0],bon[1][4]], [clip[1][0],clip[1][4]], [cam[1][0],cam[1][4]],'RS_all_5_20_plos.png', '0.5 ', '2.0 ')

	# rocPlotAll([bon[0][0],bon[0][1]], [clip[0][0],clip[0][1]], [cam[0][0],cam[0][1]], 'NP_all_5.png')