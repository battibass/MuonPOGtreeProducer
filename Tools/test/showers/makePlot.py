import argparse
import sys
import os
import json

"""
Setup argument parser
"""

parser = argparse.ArgumentParser(description="This program takes an input JSON config and extracts plots from Tag-and-Probe ROOT files. The output consists of a plot with superimposed graphs from multiple TnP files and the according fit canvases.")
parser.add_argument("inputJsonConfig", help="Path to the input JSON config file")
parser.add_argument("-f", "--fast", default=0, action="count", help="Skip fetching and saving the fit canvases for each plot")
parser.add_argument("-v", "--verbosity", default=1, help="Increase or decrease output verbosity")
args = parser.parse_args()

"""
Parse JSON file
"""

with open(args.inputJsonConfig, 'r') as f:
    config = json.loads(f.read())

"""
Go through plots defined in config JSON
"""

from ROOT import * # import this here, otherwise it overwrites the argparse stuff
gROOT.SetBatch(True) # set ROOT to batch mode, this suppresses printing canvases
gROOT.ProcessLine("gErrorIgnoreLevel = 1001;") # suppress stdout pollution of canvas.Print(...)

TH1.AddDirectory(False)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

for keyPlot in config:
    if args.verbosity==1:
        print('Processing plot config: {}'.format(keyPlot))
        print('Config comment: {}'.format(config[keyPlot]['comment']))

    # Get result plots from fit canvases (in dir fit_eff_plots)

    # Get input parameters
    inputFilenames = []
    inputPlotNames = []
    inputFolders = []
    inputLegendEntries = []
    inputLabels = []
 
    for keyInputs in sorted(config[keyPlot]['inputs']):
        inputFilenames.append(config[keyPlot]['inputs'][keyInputs]['filename'])
        inputPlotNames.append(config[keyPlot]['inputs'][keyInputs]['plot'])
        inputFolders.append(config[keyPlot]['inputs'][keyInputs]['folder'])
        inputLegendEntries.append(config[keyPlot]['inputs'][keyInputs]['legendEntry'])
        inputLabels.append(config[keyPlot]['inputs'][keyInputs]['label'])

    # Get input graphs from files
    inputHistos = {}
    for iHisto in range(len(inputFilenames)):
        inputFile = TFile.Open(inputFilenames[iHisto])
        if not inputFile:
            print('[ERROR] File not found: {}'.format(inputFilenames[iHisto]))
            sys.exit()
        inputDir = inputFile.GetDirectory(inputFolders[iHisto])
        inputName = None
        for keys in inputDir.GetListOfKeys():
            if keys.GetName() == inputPlotNames[iHisto] or inputPlotNames[iHisto] == "all" :
                inputName = keys.GetName()

                if inputName.find("=") > 0 or \
                   inputDir.Get(inputName).ClassName() == "TDirectoryFile":
                    continue

                if inputPlotNames[iHisto] == "all" :
                    nameTag = inputName
                else :
                    nameTag = "histo"

                if not inputHistos.has_key(nameTag):
                    inputHistos[nameTag] = []
            
                print('Load plot \'{}\': {}'.format(inputFolders[iHisto], inputName))
                print inputDir.Get(inputName).ClassName()

                histo = inputDir.Get(inputName).Clone(str(iHisto))
                inputHistos[nameTag].append(histo)

    # Set line color and marker style for each graph using given maps from config
    colorMap = config[keyPlot]['plot']['colorMap']
    markerMap = config[keyPlot]['plot']['markerMap']
    if args.verbosity==1:
        print('Using colormap: {}'.format(colorMap))
        print('Using markermap: {}'.format(markerMap))
    if len(colorMap)<len(inputFolders):
        print('[ERROR] The defined colormap has not enough entries for the number of defined input folders.')
        sys.exit()
    if len(markerMap)<len(inputFolders):
        print('[ERROR] The defined markermap has not enough entries for the number of defined input folders.')
        sys.exit()

    for histoName,histograms in inputHistos.iteritems():

        if histograms[iHisto].ClassName() != "TH2F" :
            for iHisto in range(len(histograms)):

                histograms[iHisto].SetLineColor(colorMap[iHisto])
                histograms[iHisto].SetLineWidth(2)
                histograms[iHisto].SetMarkerStyle(markerMap[iHisto])
                histograms[iHisto].SetMarkerColor(colorMap[iHisto])

        # Setup canvas with all elements
        canvas = TCanvas('canvas', 'canvas', 800, 800)

        pad = TPad('pad', 'pad', 0.01, 0.00, 1.00, 1.00)
 
        pad.SetGrid()
        pad.Draw()
        pad.cd()

        hasLogY      = False
        hasHistN     = False
        hasNoScale   = False
        hasLinearFit = False
        hasParabFit  = False
        
        
        if config[keyPlot]['plot'].has_key('option') :
            option = config[keyPlot]['plot']['option']

            if option.find("logY") > -1 :
                hasLogY = True
                option = option.replace("logY","")

            if option.find("histN") > -1 :
                hasHistN = True
                option = option.replace("histN","")

            if option.find("noScale") > -1 :
                hasNoScale = True
                option = option.replace("noScale","")

            if option.find("linearFit") > -1 :
                hasLinearFit = True
                option = option.replace("linearFit","")

            if option.find("parabFit") > -1 :
                hasParabFit = True
                option = option.replace("parabFit","")

        if hasLogY :
            pad.SetLogy()

        plotX = config[keyPlot]['plot']['x']
        plotY = config[keyPlot]['plot']['y']
        # Generate superimposed graph using TMultiHisto
    
        for iHisto in range(len(histograms)):

            if histograms[iHisto].ClassName() == "TH2F" :
                pad.SetLogz()
                histograms[iHisto].Draw('colz')
            else:
                if iHisto == 0 :
                    histograms[iHisto].Draw('')
                else :
                    if hasHistN :
                        histograms[iHisto].Draw('sameHIST')
                    else:
                        histograms[iHisto].Draw('same')

            canvas.Update()

            if histograms[iHisto].ClassName() == "TEfficiency" :
                histo = histograms[iHisto].GetPaintedGraph()
            else :
                histo = histograms[iHisto]

            if hasLinearFit :
               func = TF1("myLine","pol1")
               func.SetRange(plotX[0], plotX[1])
               histo.Fit("myLine","RF")
               func.SetLineColor(colorMap[iHisto])
               func.DrawCopy("same")
               par0 = float("{0:.4f}".format(func.GetParameter(0)))
               par1 = float("{0:.7f}".format(func.GetParameter(1)))
               parErr0 = float("{0:.4f}".format(func.GetParError(0)))
               parErr1 = float("{0:.7f}".format(func.GetParError(1)))
               print "[makePlot.py] Fit (line) ", histograms[iHisto].GetName(), \
                   "par0 =", par0, "parErr0 =", parErr0, \
                   "par1 =", par1, "parErr1 =", parErr1, \
               
            if hasParabFit :
               func = TF1("myParab","pol2")
               func.SetRange(plotX[0], plotX[1])
               func.SetParameter(0,0.005)
               func.SetParameter(1,0.0001)
               func.SetParameter(2,0.000000001)
               func.SetParLimits(0,0.,0.02)
               func.SetParLimits(1,0.,0.002)
               func.SetParLimits(2,-1e-7,0.)
               histo.Fit("myParab","RF")
               func.SetLineColor(colorMap[iHisto])
               func.DrawCopy("same")
               par0 = float("{0:.4f}".format(func.GetParameter(0)))
               par1 = float("{0:.7f}".format(func.GetParameter(1)))
               par2 = float("{0:.9f}".format(func.GetParameter(2)))
               parErr0 = float("{0:.4f}".format(func.GetParError(0)))
               parErr1 = float("{0:.7f}".format(func.GetParError(1)))
               parErr2 = float("{0:.9f}".format(func.GetParError(2)))
               print "[makePlot.py] Fit (parab) ", histograms[iHisto].GetName(), \
                   "par0 =", par0, "parErr0 =", parErr0, \
                   "par1 =", par1, "parErr1 =", parErr1, \
                   "par2 =", par2, "parErr2 =", parErr2




            if histograms[iHisto].ClassName() != "TH1F" :
                histo.GetXaxis().SetRangeUser(plotX[0], plotX[1])
                histo.GetYaxis().SetRangeUser(plotY[0], plotY[1])
                histo.GetXaxis().SetTitle(plotX[2])
                histo.GetYaxis().SetTitle(plotY[2])
            else:
                if not hasNoScale :
                    histo.Scale(1./histo.Integral())
                histo.GetXaxis().SetRangeUser(plotX[0], plotX[1])

                if hasLogY :
                    histo.GetYaxis().SetRangeUser(plotY[0], plotY[1])
                else :
                    histo.GetYaxis().SetRangeUser(0.0, histograms[iHisto].GetMaximum() * 1.5)

                histo.GetXaxis().SetTitle(histoName)
                histo.GetXaxis().SetTitle(plotX[2])
                histo.GetYaxis().SetTitle(plotY[2])


            histo.GetXaxis().SetLabelSize(22)
            histo.GetXaxis().SetTitleFont(63)
            histo.GetXaxis().SetLabelFont(43)
            histo.GetXaxis().SetTitleSize(22)
            histo.GetXaxis().SetLabelSize(20)
            histo.GetXaxis().SetTitleOffset(1.2)
           
            histo.GetYaxis().SetLabelSize(22)
            histo.GetYaxis().SetTitleFont(63)
            histo.GetYaxis().SetLabelFont(43)
            histo.GetYaxis().SetTitleSize(22)
            histo.GetYaxis().SetLabelSize(20)
            histo.GetYaxis().SetTitleOffset(1.5)
            
            canvas.Update()


        canvas.Update()

        canvas.cd()
        leg = TLegend(0.49, 0.67, 0.75+0.1, 0.80)

        if histograms[iHisto].ClassName() != "TH2F" :
            for iHisto in range(len(histograms)):
                leg.AddEntry(histograms[iHisto], inputLegendEntries[iHisto], 'LP')

            leg.SetBorderSize(1)
            leg.SetTextFont(43)
            leg.SetTextSize(20)
            leg.Draw()

        canvas.cd()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextFont(61)
        latex.SetTextSize(0.06)
        latex.DrawLatex(0.16, 0.82, config[keyPlot]['plot']['logo'][0])
        latex.SetTextFont(52)
        latex.SetTextSize(0.04)
        latex.SetTextAlign(11);
        latex.DrawLatex(0.16, 0.77, config[keyPlot]['plot']['logo'][1])
        latex.SetTextFont(42)
        latex.SetTextSize(0.038)
        latex.SetTextAlign(31);
        latex.DrawLatex(0.90, 0.91, config[keyPlot]['plot']['caption'])
        latex.SetTextAlign(11);
        latex.SetTextColor(1)
        latex.SetTextFont(43)
        latex.SetTextSize(20)
        latex.DrawLatex(0.49, 0.82, config[keyPlot]['plot']['legendTitle'])
        canvas.Update()
    
        # Save plot
        outputDirectory = config[keyPlot]['output']['directory']
        if args.verbosity==1:
            print('Output directory: {}'.format(outputDirectory))
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)
        for fileType in config[keyPlot]['output']['fileType']:

            if config[keyPlot]['output']['filenamePlot'] == "all":
                 fileNamePlot = histoName
            else:
                 fileNamePlot = config[keyPlot]['output']['filenamePlot']
            
            canvas.SaveAs(os.path.join(outputDirectory,fileNamePlot+'.'+fileType))

    if args.verbosity==1:
        print('')
