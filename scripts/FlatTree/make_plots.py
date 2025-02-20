import ROOT
from array import array
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

geometries = ["rmms"]#["2Coils", "4Coils", "5Coils", "SPY"]
recs = ["contained", "rmms", "none"]
nutypes = ["neutrino", "antineutrino" ]
rnames = ["LAr contained", "SSRI", "All"]
colors = [ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+2]
#topdir = "/dune/data/users/gsdavies"
topdir = "/dune/data/users/gsdavies/ssri/inc_lar/"

hSig = ROOT.TH1D( "sigma", ";#mu KE (GeV)", 20, 0., 5. )
hRMS = ROOT.TH1D( "rms", ";#mu KE (GeV)", 20, 0., 5. )

def CornerLabel( string ):
    label = ROOT.TLatex()
    label.SetTextColor(ROOT.kGray+1)
    label.SetNDC()
    label.SetTextSize(2/30.)
    label.SetTextAlign(11)
    label.DrawLatex(.1, .93, string)

def setNormalColors():
    stops = array( 'd', [ 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000 ] )
    red   = array( 'd', [0./255.,  61./255.,  89./255., 122./255., 143./255., 160./255., 185./255., 204./255., 231./255.] )
    green = array( 'd', [0./255.,   0./255.,   0./255.,   0./255.,  14./255.,  37./255.,  72./255., 132./255., 235./255.] )
    blue  = array( 'd', [0./255., 140./255., 224./255., 144./255.,   4./255.,   5./255.,   6./255.,   9./255.,  13./255.] )
    ROOT.TColor.CreateGradientColorTable(9, stops, red, green, blue, 255, 1);

def setRedBlue():
    stops = array( 'd', [0., 0.5, 1.] )
    red   = array( 'd', [0., 1., 1.] )
    green = array( 'd', [0., 1., 0.] )
    blue  = array( 'd', [1., 1., 0.] )
    ROOT.TColor.CreateGradientColorTable(3, stops, red, green, blue, 999, 1);

def plotVXY( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")
    c.Print( "testplots/vxy_%s_test.png" % (geom) )

def plotHitPos( ohdenom, geom, xyz, birth ):
    hdenom = ohdenom.Clone( "clone" )

    c = ROOT.TCanvas()
    hdenom.Draw("colz")
    if birth == True:
        c.Print( "testplots/birth%s_%s.png" % (xyz, geom) )
    else:
        c.Print( "testplots/death%s_%s.png" % (xyz, geom) )

def plotED( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")

    #prof = hdenom.ProfileX();
    #prof.SetLineColor(ROOT.kRed)
    #prof.Fit("pol1","","",0,1);
    #fitter = prof.FindObject("pol1")
    #fitter.SetLineColor(ROOT.kYellow+2)
    #hdenom.Draw("colz same")
    #c.SetLogx()

    hdenom.GetXaxis().CenterTitle()
    hdenom.GetYaxis().CenterTitle()
    hdenom.GetYaxis().SetTitle("Reconstructed Track Length (g/cm^{2})")
    #hdenom.GetXaxis().SetTitle("Muon Momentum (GeV) [SSRI]")
    hdenom.GetXaxis().SetTitle("Muon Kinetic Energy (GeV) [SSRI]")

    #c.RedrawAxis()
    #hdenom.Draw("colz same")
    c.Print( "testplots/MOM-EDP_%s_test.png" % (geom) )

def plotED2( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")

    #prof = hdenom.ProfileX();
    #prof.SetLineColor(ROOT.kRed)
    #prof.Fit("pol1","","",0,5);
    #fitter = prof.FindObject("pol1")
    #fitter.SetLineColor(ROOT.kYellow+2)
    #hdenom.Draw("colz same")
    #c.SetLogx()

    hdenom.GetXaxis().CenterTitle()
    hdenom.GetYaxis().CenterTitle()
    hdenom.GetYaxis().SetTitle("Reconstructed Track Length (g/cm^{2})")
    hdenom.GetXaxis().SetTitle("Muon Momentum (GeV) [SSRI]")

    #c.SetLogx()
    c.Print( "testplots/MOMorig-EDP_%s_test.png" % (geom) )

def plotEDTrk( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")
    c.SetLogx()
    c.Print( "testplots/edep-trk_%s_test.png" % (geom) )

def plotMomRes( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c1 = ROOT.TCanvas()
    hdenom.Draw("colz")
    #c1.Print( "testplots/mom-recotrue_%s_test.png" % (geom) )
    
    c = ROOT.TCanvas()
    for b in range(1,21):
        h = hdenom.ProjectionY("py_%d"%b, b, b)

        hRMS.SetBinContent( b, h.GetRMS() )
        maxP = h.GetBinCenter(h.GetMaximumBin() )
        print("maxP", maxP)
        if h.GetSumOfWeights() != 0:
            h.Fit( "gaus", "Q", "Q", maxP-0.2, maxP+0.2 )
            print("sigma ", h.GetFunction("gaus").GetParameter(2))
            hSig.SetBinContent( b, h.GetFunction("gaus").GetParameter(2) )
            h.Draw("hist")
            h.GetFunction("gaus").SetLineColor(2)
            h.GetFunction("gaus").Draw("same")
            #c.Print( "testplots/residual_%1.1f_%1.1f.png" % (hdenom.GetXaxis().GetBinLowEdge(b), hdenom.GetXaxis().GetBinLowEdge(b+1)) )
    
    hRMS.SetMinimum(0.)
    hRMS.SetLineColor(2)
    hRMS.Draw("hist")
    hSig.Draw("hist same")
    c.Print( "testplots/rms_sigma.png" )
    # hdenom.SetLineColor(ROOT.kAzure+2)
    # gaussFit = ROOT.TF1("gaussfit","gaus",-1.,1.)
    # hdenom.Fit(gaussFit,"E")

    # latex = ROOT.TLatex()
    # latex.SetNDC()
    # latex.SetTextSize(0.03)
    # mean = gaussFit.GetParameter(1)
    # sigma = gaussFit.GetParameter(2)
    # ndof = gaussFit.GetNDF()
    # chi2 = gaussFit.GetChisquare()
    
    # rms = hdenom.GetRMS()
    # print "RMS: ",rms

    # latex.DrawText(0.6,0.7, "Mean = %.1f"%(mean))
    # latex.DrawText(0.6,0.65, "Sigma = %.1f"%(sigma))
    # latex.DrawText(0.6,0.6, "#chi^{2}/ndof = %.1f/%d = %.1f"%(chi2,ndof,chi2/ndof))
    # latex.DrawText(0.6,0.55, "RMS = %.1f"%(rms))
    
    #fitter = .FindObject("")
    #fitter.SetLineColor(ROOT.kBlack)

    #hdenom.Draw("e")
    #prof = hdenom.ProfileX()
    #prof.SetLineColor(ROOT.kRed+2)
    #prof.Draw()

def plotMomRes2( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")

    #prof = hdenom.ProfileX("name2",-4,4,"s")
    #prof.SetLineColor(ROOT.kRed)
    
    #prof.Fit("pol1","","",0,5);
    #fitter = prof.FindObject("pol1")
    #fitter.SetLineColor(ROOT.kYellow+2)
    c.Update()
    c.Print( "testplots/mom-recotrue_v_true_%s_test.png" % (geom) )
    c.Print( "testplots/mom-recotrue_v_true_%s_test.root" % (geom) )

def plotAngle( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")
    c.Print( "testplots/angle_%s_test.png" % (geom) )

def plotLArFactor( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("hist")
    c.Print( "testplots/larfactor_%s_test.png" % (geom) )

def plotScint( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("hist")
    c.Print( "testplots/scint_%s_test.png" % (geom) )

def plotVYZ( ohdenom, geom ):
    hdenom = ohdenom.Clone( "clone" )
    
    c = ROOT.TCanvas()
    hdenom.Draw("colz")
    c.Print( "testplots/vyz_%s_test.png" % (geom) )

def plotMom( ac, ssri, con, geom ):
    hac = ac.Clone( "clone1" )
    hssri = ssri.Clone( "clone2" )
    hcon = con.Clone( "clone3" )
    
    c = ROOT.TCanvas()
    hac.Draw("hist")
    hac.SetLineColor(ROOT.kRed+2)
    hssri.Draw("hist same")
    hssri.SetLineColor(ROOT.kAzure+2)
    hcon.Draw("hist same")
    hcon.SetLineColor(ROOT.kGreen+2)

    hac.GetXaxis().CenterTitle()
    hac.GetYaxis().CenterTitle()
    hac.GetYaxis().SetTitle("No. of events")
    hac.GetXaxis().SetTitle("Muon Momentum (GeV)")

    leg = ROOT.TLegend( 0.2, 0.75, 0.8, 0.846 )
    leg.SetNColumns(3)
    leg.AddEntry( hac, rnames[2], "f" )
    leg.AddEntry( hssri, rnames[1], "f" )
    leg.AddEntry( hcon, rnames[0], "f" )
    leg.Draw()
    c.Print( "testplots/momentum_comp_%s_test.png" % (geom) )


def plotAcc2D( ohNum, ohDenom, geom, name ):
    hNum = [ohNum[i].Clone("num_%s" % r) for i,r in enumerate(recs)]
    hDenom = ohDenom.Clone( "clone" )

    hNum[1].Add( hNum[0] )
    hNum[1].Add( hNum[2] )

    for i in range(3):
        hNum[i].Divide(hDenom)
        hNum[i].SetMinimum(0.)
        hNum[i].SetMaximum(1.)

    c = ROOT.TCanvas()
    hNum[0].Draw("colz")
    c.Print( "testplots/acc2d_%s_%s_contained.png" % (name,geom) )
    c.Print( "testplots/acc2d_%s_%s_contained.pdf" % (name,geom) )

    hNum[1].Draw("colz")
    c.Print( "testplots/acc2d_%s_%s_contained+rmms.png" % (name,geom) )
    c.Print( "testplots/acc2d_%s_%s_contained+rmms.pdf" % (name,geom) )

    hNum[2].Draw("colz")
    c.Print( "testplots/acc2d_%s_%s_contained+rmms+none.png" % (name,geom) )
    c.Print( "testplots/acc2d_%s_%s_contained+rmms+none.pdf" % (name,geom) )

def plotELen2D( ohNum, ohDenom, geom, name ):
    hNum = [ohNum[i].Clone("numlen_%s" % r) for i,r in enumerate(recs)]
    hDenom = ohDenom.Clone( "clone" )

    c = ROOT.TCanvas()
    hNum[0].Draw("colz")
    c.Print( "testplots/elen2d_%s_%s_contained.png" % (name,geom) )
    c.Print( "testplots/elen2d_%s_%s_contained.pdf" % (name,geom) )

    hNum[1].Draw("colz")
    c.Print( "testplots/elen2d_%s_%s_contained+rmms.png" % (name,geom) )
    c.Print( "testplots/elen2d_%s_%s_contained+rmms.pdf" % (name,geom) )

    hNum[2].Draw("colz")
    c.Print( "testplots/acc2d_%s_%s_none.png" % (name,geom) )
    c.Print( "testplots/acc2d_%s_%s_none.pdf" % (name,geom) )

def plotAccRatio( ohNum, ohDenom, ohNum2, ohDenom2, geom ):
    hNum = ohNum[0].Clone("rnum")
    hNum.Add( ohNum[1] )
    #hNum.Add( ohNum[2] )
    hDenom = ohDenom.Clone( "rclone" )
    hNum1 = ohNum2[0].Clone("rnum2")
    hNum1.Add( ohNum2[1] )
    #hNum2.Add( ohNum2[2] )
    hDenom1 = ohDenom2.Clone( "r2clone" )

    hNum.Divide( hDenom )
    hNum1.Divide( hDenom1 )
    hNum.Divide( hNum1 )

    c = ROOT.TCanvas()
    hNum.SetMinimum(0.7)
    hNum.SetMaximum(1.3)
    hNum.Draw("colz")
    c.Print( "testplots/accRatio_%s.png" % geom )
    c.Print( "testplots/accRatio_%s.pdf" % geom )

def plotAccRatioAll( ohNum, ohDenom, idx ):
    c = ROOT.TCanvas()
    cccolors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kGreen+2]
    hRatios = [None for g in geometries]
    leg = ROOT.TLegend( 0.2, 0.16, 0.8, 0.22 )
    leg.SetNColumns(4)
    for g,geom in enumerate(geometries):
        hRatios[g] = ohNum[g][0].ProjectionX("aaa%s" % geom, 1, 2)
        hRatios[g].Add( ohNum[g][1].ProjectionX("aaaa%s" % geom, 1, 2))
        #hRatios[g].Add( ohNum[g][2].ProjectionX("aaaaa%s" % geom, 1, 2))
        hD1d = ohDenom[g].ProjectionX("bbb%s" % geom, 1, 2)
        hN1d2 = ohNum[idx][0].ProjectionX("ccc%s" % geom, 1, 2)
        hN1d2.Add( ohNum[idx][1].ProjectionX("cccc%s" % geom, 1, 2))
        #hN1d2.Add( ohNum[idx][2].ProjectionX("ccccc%s" % geom, 1, 2))
        hD1d2 = ohDenom[idx].ProjectionX("ddd%s" % geom, 1, 2)

        hRatios[g].Divide( hD1d )
        hN1d2.Divide( hD1d2 )
        hRatios[g].Divide( hN1d2 )
        hRatios[g].SetLineColor(cccolors[g])
        hRatios[g].Rebin(2)
        hRatios[g].Scale(0.5)
        leg.AddEntry( hRatios[g], geom, "l" )

    hRatios[0].SetMinimum(0.9)
    hRatios[0].SetMaximum(1.05)
    hRatios[0].SetTitle( ";Muon kinetic energy (GeV);Acc. ratio to %s" % geometries[idx] )
    hRatios[0].Draw("hist")
    for i in range(1,4):
        hRatios[i].Draw("hist same")
    leg.Draw()
    c.Print( "testplots/accRatioAll.png" )
    c.Print( "testplots/accRatioAll.pdf" )

def plotSignDist( ohNum, geom ):
    hNum = [ohNum[t].Clone("num_%s" % nutype) for t,nutype in enumerate(nutypes)]
    
    c = ROOT.TCanvas()
    hNum[0].SetLineColor(ROOT.kAzure+2)
    hNum[1].SetLineColor(ROOT.kRed+2)
    hNum[0].Draw("hist")
    integral = hNum[0].Integral(50,100)
    totintegral = hNum[0].Integral(-1,-1)
    print("INTEGRAL: ", integral)
    print("tot INTEGRAL: ", totintegral)
    hNum[1].Draw("hist same")
    leg = ROOT.TLegend( 0.2, 0.75, 0.8, 0.846 )
    leg.SetNColumns(2)

    leg.AddEntry( hNum[0], "Neutrino", "l" )
    leg.AddEntry( hNum[1], "Antineutrino", "l" )


    leg.Draw()
    CornerLabel("Neutrino Beam Mode")
    c.Print("testplots/signdist_%s.png" % geom )

def plotSignDistMom( ohNum, geom ):
    hNum = [ohNum[t].Clone("num_%s" % nutype) for t,nutype in enumerate(nutypes)]
    
    c = ROOT.TCanvas()
    hNum[0].Draw("colz")
    CornerLabel("Neutrino Beam Mode")
    c.Print("testplots/signdistmom_%s_neutrino.png" % geom )

    hNum[1].Draw("colz")
    CornerLabel("Neutrino Beam Mode")
    c.Print("testplots/signdistmom_%s_antineutrino.png" % geom )
    #leg = ROOT.TLegend( 0.2, 0.75, 0.8, 0.846 )
    #leg.SetNColumns(2)

    #leg.AddEntry( hNum[0], "Neutrino", "l" )
    #leg.AddEntry( hNum[1], "Antineutrino", "l" )


def plotAcc1D( ohNum, ohDenom, geom ):
    hNum = [ohNum[i].ProjectionX("num_%s" % r, 1, 2) for i,r in enumerate(recs)]
    hDenom = ohDenom.ProjectionX("proj", 1, 2)

    stk = ROOT.THStack( "stk", "#theta_{#mu} < 20 degrees;Muon kinetic energy (GeV);Acceptance" )
    leg = ROOT.TLegend( 0.2, 0.75, 0.8, 0.846 )
    leg.SetNColumns(3)

    for i in [0, 1, 2]:
        hNum[i].Divide(hDenom)
        hNum[i].SetLineColor(colors[i])
        hNum[i].SetFillColor(colors[i])
        hNum[i].SetFillStyle(1)
        hNum[i].SetMinimum(0.)
        hNum[i].SetMaximum(1.2)
        stk.Add( hNum[i] )
        leg.AddEntry( hNum[i], rnames[i], "f" )

    c = ROOT.TCanvas()
    stk.Draw("hist")
    stk.SetMaximum(1.1)
    stk.Draw("hist")
    leg.Draw()
    c.Print( "testplots/acc1d_KE_stk_%s.png" % geom )
    c.Print( "testplots/acc1d_KE_stk_%s.pdf" % geom )

    # non-stack
    leg = ROOT.TLegend(0.6, 0.2, 0.846, 0.5)
    hNumAll = hNum[0].Clone( "numall" )
    hNumAll.SetLineColor(1)
    hNumAll.SetFillStyle(0)
    hNumAll.Add( hNum[1] )
    hNumAll.Add( hNum[2] )
    leg.AddEntry( hNumAll, "Total", "l" )
    hNum[0].SetMinimum(0.)
    hNum[0].SetMaximum(1.)
    hNum[0].SetTitle( "#theta_{#mu} < 20 degrees;Muon kinetic energy (GeV);Acceptance" )
    for i in [0, 1, 2]:
        hNum[i].SetFillStyle(0)
        leg.AddEntry( hNum[i], rnames[i], "l" )
        hNum[i].Draw("hist" if not i else "hist same")
    hNumAll.Draw("hist same")
    leg.Draw()
    c.Print( "testplots/acc1d_KE_%s.png" % geom )
    c.Print( "testplots/acc1d_KE_%s.pdf" % geom )

def plotX( hNum, hDenom, geom ):
    c = ROOT.TCanvas()

    energies = [0.7, 0.9, 1.1, 1.3]
    hRatios = [None for e in energies]
    leg = ROOT.TLegend( 0.2, 0.16, 0.8, 0.22 )
    leg.SetNColumns(len(energies))
    for e,energy in enumerate(energies):
        hRatios[e] = hNum.ProjectionY("asdf_%d" % e, hNum.GetXaxis().FindBin(energy), hNum.GetXaxis().FindBin(energy))
        hRatios[e].Divide( hDenom.ProjectionY("asdfasdf_%d" % e, hNum.GetXaxis().FindBin(energy), hNum.GetXaxis().FindBin(energy)) )
        hRatios[e].Rebin(4)
        hRatios[e].Scale(0.25)
        hRatios[e].SetLineColor(colors[e])
        leg.AddEntry( hRatios[e], "T_{#mu} = %1.1f" % energy, "l" )
        if not e:
            hRatios[e].SetTitle( ";X position (cm);#theta_{#mu} < 30^{#circ} acceptance" )
            hRatios[e].SetMinimum(0.3)
            hRatios[e].SetMaximum(0.8)
            hRatios[e].Draw("hist")
        else:
            hRatios[e].Draw("hist same")
    leg.Draw()
    c.Print( "testplots/accX_%s.png" % geom )
    c.Print( "testplots/accX_%s.pdf" % geom )

    hRatio = hNum.ProjectionY("easdf", 4, 8)
    hRatio.Divide( hDenom.ProjectionY("asdfasdf", 4, 8) )
    hRatio.Draw("hist")
    c.Print( "testplots/accX_dip_%s.png" % geom )
    c.Print( "testplots/accX_dip_%s.pdf" % geom )

    hNum.Divide(hDenom)
    hNum.SetMinimum(0.3)
    hNum.SetMaximum(0.8)
    hNum.Draw("colz")
    c.Print( "testplots/accX2D_%s.png" % geom )
    c.Print( "testplots/accX2D_%s.pdf" % geom )

if __name__ == "__main__":

    tf = ROOT.TFile( "%s/dumpSSRI-v5-histos.root" % topdir )
    #tf = ROOT.TFile( "%s/filled_fifty.root" % topdir )

    hDenom = [None for g in geometries]
    hNum = [[None for r in recs] for g in geometries]
    hSignDist = [[None for t in nutypes] for g in geometries]
    hSignDistMom = [[None for t in nutypes] for g in geometries]
    hNumLen = [[None for r in recs] for g in geometries]
    hDenomQ2 = [None for g in geometries]
    hNumQ2 = [[None for r in recs] for g in geometries]
    hDenomX = [None for g in geometries]
    hNumX = [None for g in geometries]
    hVXY = [None for g in geometries]
    hVYZ = [None for g in geometries]
    hMuonBirthXZ = [None for g in geometries]
    hMuonBirthYZ = [None for g in geometries]
    hMuonDeathXZ = [None for g in geometries]
    hMuonDeathYZ = [None for g in geometries]
    hED = [None for g in geometries]
    hED2 = [None for g in geometries]
    hEDTrk = [None for g in geometries]
    hScint = [None for g in geometries]
    hMomSSRI = [None for g in geometries]
    hMomAC = [None for g in geometries]
    hMomCON = [None for g in geometries]
    hMomRes = [None for g in geometries]
    hMomRes2 = [None for g in geometries]
    hAngle = [None for g in geometries]
    hLArFactor = [None for g in geometries]

    setNormalColors()
    for g,geom in enumerate(geometries):
        hDenom[g] = tf.Get( "denom_%s" % geom )
        hDenomQ2[g] = tf.Get( "denomQ2_%s" % geom )
        hDenomX[g] = tf.Get( "denomX_%s" % geom )
        hNumX[g] = tf.Get( "numX_%s" % geom )
        hVXY[g] = tf.Get( "denomvxy_%s" % geom )
        hVYZ[g] = tf.Get( "denomvyz_%s" % geom )
        hMuonBirthXZ[g] = tf.Get( "birthxz_%s" % geom )
        hMuonBirthYZ[g] = tf.Get( "birthyz_%s" % geom )
        hMuonDeathXZ[g] = tf.Get( "deathxz_%s" % geom )
        hMuonDeathYZ[g] = tf.Get( "deathyz_%s" % geom )
        hED[g] = tf.Get( "MOM-EDEP_%s" % geom )
        hED2[g] = tf.Get( "MOMorig-EDEP_%s" % geom )
        hEDTrk[g] = tf.Get( "edep-trk_%s" % geom )
        hScint[g] = tf.Get( "scint_%s" % geom )
        hMomAC[g] = tf.Get( "momAC_%s" % geom )
        hMomSSRI[g] = tf.Get( "momSSRI_%s" % geom )
        hMomCON[g] = tf.Get( "momCON_%s" % geom )
        hMomRes[g] = tf.Get( "momRES_%s" % geom )
        hMomRes2[g] = tf.Get( "momRES2_%s" % geom )
        hAngle[g] = tf.Get( "angle_%s" % geom )
        hLArFactor[g] = tf.Get( "larfactor_%s" % geom )
        for r,rec in enumerate(recs):
            hNum[g][r] = tf.Get( "num_%s_%s" % (geom, rec) )
            hNumLen[g][r] = tf.Get( "numlen_%s_%s" % (geom, rec) )
            hNumQ2[g][r] = tf.Get( "numQ2_%s_%s" % (geom, rec) )

        for t, nutype in enumerate(nutypes):
            hSignDist[g][t] = tf.Get("signdist_%s_%s" % (geom, nutype) )
            hSignDistMom[g][t] = tf.Get("signdistmom_%s_%s" % (geom, nutype) )

        plotSignDist( hSignDist[g], geom)
        plotSignDistMom( hSignDistMom[g], geom)
            
        plotAcc2D( hNum[g], hDenom[g], geom, "etheta" )
        plotELen2D( hNumLen[g], hDenom[g], geom, "elen" )
        plotAcc2D( hNumQ2[g], hDenomQ2[g], geom, "Q2" )
        plotAcc1D( hNum[g], hDenom[g], geom )
        plotX( hNumX[g], hDenomX[g], geom )
        plotVXY( hVXY[g], geom )
        plotVYZ( hVYZ[g], geom )
        plotHitPos( hMuonBirthXZ[g], geom, "XZ", True )
        plotHitPos( hMuonBirthYZ[g], geom, "YZ", True )
        plotHitPos( hMuonDeathXZ[g], geom, "XZ", False )
        plotHitPos( hMuonDeathYZ[g], geom, "YZ", False )
        plotED( hED[g], geom )
        plotED2( hED2[g], geom )
        plotEDTrk( hEDTrk[g], geom )
        plotMomRes( hMomRes[g], geom )
        plotMomRes2( hMomRes2[g], geom )
        plotScint( hScint[g], geom )
        plotMom( hMomAC[g], hMomSSRI[g], hMomCON[g], geom )
        plotAngle( hAngle[g], geom)
        plotLArFactor( hLArFactor[g], geom)

        #setRedBlue()
    #for g,geom in enumerate(geometries):
    #    plotAccRatio( hNum[g], hDenom[g], hNum[2], hDenom[2], geom )

    #plotAccRatioAll( hNum, hDenom, 2 )

