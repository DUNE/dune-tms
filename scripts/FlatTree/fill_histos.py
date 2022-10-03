import ROOT
from math import cos, atan2, sqrt

geometries = ["rmms"]#["2Coils", "4Coils", "5Coils", "SPY"]
#topdir = "/dune/data/users/gsdavies/"
#topdir="/dune/data/users/gsdavies/ssri/inc_lar/"
topdir="/dune/app/users/cwret/SSRI/dune-ssri/ana"
recs = ["contained", "rmms", "none"]
nutype = ["neutrino", "antineutrino" ]

def loop( tree, hNum, hNumLen, hDenom, hNumX, hDenomX, hNumQ2, hDenomQ2, hVXY, hVYZ, hMuonBirthXZ, hMuonBirthYZ, hMuonDeathXZ, hMuonDeathYZ, hAngle, hED, hED2, hEDTrk, hScint, hMomAC, hMomSSRI, hMomCON, hMomRes, hMomRes2, hLArFactor, hSignDist, hSignDistMom ):

    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("vtx", 1)
    tree.SetBranchStatus("ievt", 1)
    tree.SetBranchStatus("Ev", 1)
    tree.SetBranchStatus("lepPdg", 1)
    tree.SetBranchStatus("lepKE", 1)
    tree.SetBranchStatus("p3lep", 1)
    tree.SetBranchStatus("lepDeath", 1)
    tree.SetBranchStatus("muonDeath", 1)
    tree.SetBranchStatus("muonBirth", 1)
    tree.SetBranchStatus("muScintLen", 1)
    tree.SetBranchStatus("muLArLen", 1)
    tree.SetBranchStatus("muScintEnergy", 1)
    tree.SetBranchStatus("muonReco", 1)
    tree.SetBranchStatus("muScintLenOld", 1)
    tree.SetBranchStatus("muonExitPt", 1)
    tree.SetBranchStatus("muonExitMom", 1)
    tree.SetBranchStatus("lepRMMS", 1)
    tree.SetBranchStatus("muonExitKE", 1)
    tree.SetBranchStatus("rmmsKE", 1)

    beam_angle = ROOT.TVector3( 0., 0., 1. )
    beam_angle.RotateX( 0.101 )

    N = tree.GetEntries()
    n = 0
    total_selected = 0.
    total = 0.
    print("entries: ", N)
    for entry in tree:
        if n % 1000 == 0:
            print("  event %d of %d..." % (n,N))
        n += 1
        print("event ", n)

        # True numu CC (or rather just mu+/-)
        if abs(entry.lepPdg) != 13: continue
        print("passed lepton req")

        # vertex cut (LAr fiducial)
        if abs(entry.vtx[0]) > 300. or abs(entry.vtx[1]) > 100. or entry.vtx[2] < 50. or entry.vtx[2] > 350.: continue

        # true muon momentum
        mu = ROOT.TVector3( entry.p3lep[0], entry.p3lep[1], entry.p3lep[2] )
        # Muon momentum when exiting active LAr volume
        muExit = ROOT.TVector3( entry.muonExitMom[0], entry.muonExitMom[1], entry.muonExitMom[2] )
        theta = mu.Angle(beam_angle)
        thetadeg = (180./3.1416)*theta

        # Lepton Kinetic energy (convert to MeV w/ adjustment for muon mass)
        Elep = entry.lepKE/1000. + 0.105658
        plep = mu.Mag()/1000.
        plepexit = muExit.Mag()/1000.
        Q2 = 2*entry.Ev*(Elep - plep*cos(theta)) - 0.105658**2

        
        hDenom.Fill( entry.lepKE/1000., thetadeg )
        hDenomQ2.Fill( entry.Ev, Q2 )

        #hVXY.Fill( entry.vtx[0], entry.vtx[1] )
        #hVYZ.Fill( entry.vtx[2], entry.vtx[1] )
       
        hMomAC.Fill(plep)
        larKE = entry.lepKE - entry.muonExitKE
        rec = False
        # LAr contained
        if entry.muonReco == 1: 
            hNum[0].Fill( entry.lepKE/1000., thetadeg )
            hNumLen[0].Fill( entry.lepKE/1000., entry.muScintLen )
            hNumQ2[0].Fill( entry.Ev, Q2 )
            hMomCON.Fill( plep )
            rec = True
        # endpoint not in active material (i.e. uncontained)
        elif entry.muonReco == 0:
            hNum[2].Fill( entry.lepKE/1000., thetadeg )
            hNumQ2[2].Fill( entry.Ev, Q2 )
            rec = False
        # SSRI contained
        elif entry.muonReco == 2:
            #hNum[1].Fill( entry.lepKE/1000., thetadeg )
            #hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
            #hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
            #hNumQ2[1].Fill( entry.Ev, Q2 )
            rec = True
            # x -340
        match = False
        signed_dist = 0.

        # TODO: Make these easier to understand
        # muonBirth = first muon hit in SSRI
        # muonDeath = last muon hit in SSRI
        # Values from the geometry
        # Here we're cutting so we only have front-entering (Z < 733. cm)
        # i.e. nothing from the side, only the first plane
        # SSRI active cuts:
        #    Endpoint Z < 1300 cm
        #    -260 < Y < -25 cm
        #    -300 < X < 300 cm and non-zero KE
        if entry.muonBirth[2] < 733. and entry.muonDeath[2] < 1300. and entry.muonDeath[1] < 25. and entry.muonDeath[1] > -260. and entry.muonDeath[0] < 300. and entry.muonDeath[0] > -300. and entry.muonExitKE > 0.:
            #if entry.muonReco == 2 and entry.lepPdg == -13:
            #if (entry.muonDeath[0] < 165. and entry.muonDeath[0] > 10. or entry.muonDeath[0] > -165. and entry.muonDeath[0] < -10.):# and (entry.muonBirth[0] < 170. and entry.muonBirth[0] > 10. or entry.muonBirth > -170. and entry.muonBirth < -10.):
            # Now we're refining the selection further to only select
            # the stips where the uniform field exists (avoiding the
            # gappy regions in XZ space
            # (hmm, should also probably restrict in XY?
            # muon endpoint:  -165 < X < -10, 10 < X < 165 cm
            if (entry.muonDeath[0] < 165. and entry.muonDeath[0] > 10. or entry.muonDeath[0] > -165. and entry.muonDeath[0] < -10.):# and (entry.muonBirth[0] < 170. and entry.muonBirth[0] > 10. or entry.muonBirth > -170. and entry.muonBirth < -10.):

                hNum[1].Fill( entry.lepKE/1000., thetadeg )
                #hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
                hNumLen[1].Fill( entry.lepKE/1000., entry.muScintLen )
                hNumQ2[1].Fill( entry.Ev, Q2 )

                # These are used for calculating the signed distance metric
                x1 = entry.muonExitPt[0]
                z1 = entry.muonExitPt[2]
                x2 = entry.muonBirth[0]
                z2 = entry.muonBirth[2]
                x3 = entry.muonDeath[0]
                z3 = entry.muonDeath[2]
            
                num = (-(z2-z1)*x3 + (x2-x1)*z3 + (x1*z2 - z1*x2))
                denom = sqrt(((x2-x1)**2) + ((z2-z1)**2))
                # skip values where denom = 0
                if denom == 0: signed_dist = 0.
                else: signed_dist = num/denom
                #print ""
                #print " ", entry.lepPdg
                #print "Signed distance: ", signed_dist
                #Save histograms for (anti-)neutrino separately
                if entry.lepPdg == 13:
                    hSignDist[0].Fill(signed_dist)
                    hSignDistMom[0].Fill(signed_dist, plepexit)
                elif entry.lepPdg == -13:
                    hSignDist[1].Fill(signed_dist)
                    hSignDistMom[1].Fill(signed_dist, plepexit)


                total += 1
                if (entry.lepPdg == 13 and signed_dist < 0) or (entry.lepPdg == -13 and signed_dist > 0):
                    total_selected += 1
                
                    #print " Event: ", entry.ievt
                #print " Event: ", entry.lepKE

                #hLArFactor.Fill((entry.muonExitKE)/entry.muScintLen)
                hLArFactor.Fill(entry.muonExitKE/(entry.muScintLen))
                hED.Fill( entry.muonExitKE/1000., entry.muScintLen )
                hVXY.Fill( entry.lepDeath[0], entry.lepDeath[1] )
                hVYZ.Fill( entry.lepDeath[2], entry.lepDeath[1] )
                #if entry.rmmsKE/1000. > 2. and entry.muScintLen < 500.:
                hMuonDeathYZ.Fill( entry.muonDeath[2], entry.muonDeath[1] )
                hMuonDeathXZ.Fill( entry.muonDeath[2], entry.muonDeath[0] )
                hMuonBirthYZ.Fill( entry.muonBirth[2], entry.muonBirth[1] )
                hMuonBirthXZ.Fill( entry.muonBirth[2], entry.muonBirth[0] )
            
                # exit KE vs. track length
                hED.Fill( entry.muonExitKE/1000., entry.muScintLen )

                #reconstructed track length vs. "true" momentum 
                hED2.Fill( plepexit, entry.muScintLen )
                # leaving here, but this was for some testing using linear
                # energy estimation and 2D fit
                hMomRes2.Fill(plepexit, ((entry.muScintLen+5.42)/4963.)-plepexit)
                # 2D momentum resolution, reco KE vs. (reco-true/true)
                # Note: 2 MeV/cm energy loss for MIP
                # 1.9 is from approximation from LArFactor plot
                hMomRes.Fill( entry.muonExitKE/1000., (1.9*entry.muScintLen-entry.muonExitKE)/entry.muonExitKE )

                dedx = 0.
                if entry.muScintLen == 0.:
                    dedx = 0.
                else: dedx = entry.muScintEnergy/entry.muScintLen
                hEDTrk.Fill(dedx, plepexit )
                hScint.Fill(entry.muScintEnergy)
                #hAngle.Fill( entry.muScintLen, diff )
                hMomSSRI.Fill(plep)
                rec = True

        #print "TOTAL: " , total
        #print "TOTAL SEL: ", total_selected
        if thetadeg < 30.:
            hDenomX.Fill( entry.lepKE/1000., entry.vtx[0] )
            if rec:
                hNumX.Fill( entry.lepKE/1000., entry.vtx[0] )

    print(" total_selected: ", total_selected)
    print("total: ", total)

    purity = total_selected / total
    print("Purity: ", purity)

if __name__ == "__main__":

    hDenom = [None for g in geometries]

    hNum = [[None for r in recs] for g in geometries]
    hNumLen = [[None for r in recs] for g in geometries]
    hDenomQ2 = [None for g in geometries]
    hNumQ2 = [[None for r in recs] for g in geometries]
    hSignDist = [[None for t in nutype] for g in geometries]
    hSignDistMom = [[None for t in nutype] for g in geometries]
    hDenomX = [None for g in geometries]
    hNumX = [None for g in geometries]
    hVXY = [None for g in geometries]
    hVYZ = [None for g in geometries]
    #hVXZ = [None for g in geometries]
    hMuonBirthXZ = [None for g in geometries]
    hMuonBirthYZ = [None for g in geometries]
    hMuonDeathXZ = [None for g in geometries]
    hMuonDeathYZ = [None for g in geometries]
    hAngle = [None for g in geometries]
    hLArFactor = [None for g in geometries]
    hED = [None for g in geometries]
    hED2 = [None for g in geometries]
    hEDTrk = [None for g in geometries]
    hScint = [None for g in geometries]
    hMomAC = [None for g in geometries]
    hMomSSRI = [None for g in geometries]
    hMomCON = [None for g in geometries]
    hMomRes = [None for g in geometries]
    hMomRes2 = [None for g in geometries]
    for g,geom in enumerate(geometries):
        tree = ROOT.TChain("tree","tree")
        #tree.Add( "%s/MPD_%s_LAr/*.root" % (topdir, geom) )
        #tree.Add("%s/test10-60.root" % topdir )
        #tree.Add("%s/dumpSSRI-v5.root" % topdir )
        #tree.Add("%s/rhc_dump10v1.root" % topdir )
        #tree.Add("%s/test.root" % topdir )
        #tree.Add("%s/dump_cmfix_allnu.root" % topdir )
        tree.Add("%s/test_ssri_10July.root" % topdir)

        hAngle[g] = ROOT.TH2D( "angle_%s" % geom, ";Neutrino energy (GeV);Angle (RMMS) [degrees]", 200, 0., 2000., 45, -90., 90. )
        hED[g] = ROOT.TH2D( "MOM-EDEP_%s" % geom, ";Muon Kinetic Energy (GeV) [SSRI];Muon Track Length (g/cm^2)", 35, 0., 7., 40, 0., 4000. )
        hED2[g] = ROOT.TH2D( "MOMorig-EDEP_%s" % geom, ";Muon Momentum (GeV [SSRI] ;Reconstructed Muon Momentum (GeV) [SSRI]", 60, 0., 6., 100, 0., 2000. )
        hEDTrk[g] = ROOT.TH2D( "edep-trk_%s" % geom, ";Muon Momentum (GeV) [SSRI];#mu dE/dx (MeV/g/cm^{2})", 600, 0.05, 2., 500, 0.0, 5. )
        hScint[g] = ROOT.TH1D( "scint_%s" % geom, ";Energy Deposited", 100, 0., 500. )
        #hLArFactor[g] = ROOT.TH1D( "larfactor_%s" % geom, ";True KE/Track Length", 100, 1., 3. )
        hLArFactor[g] = ROOT.TH1D( "larfactor_%s" % geom, ";True KE/Track Length", 100, 0., 3. )
        hMomAC[g] = ROOT.TH1D( "momAC_%s" % geom, ";Muon Momentum (All)", 200, 0., 20. )
        hMomSSRI[g] = ROOT.TH1D( "momSSRI_%s" % geom, ";Muon Momentum (SSRI)", 200, 0., 20. )
        hMomCON[g] = ROOT.TH1D( "momCON_%s" % geom, ";Muon Momentum (ArgonCube)", 200, 0., 20. )
        hMomRes2[g] = ROOT.TH2D( "momRES2_%s" % geom, ";True Muon Momentum (GeV) [SSRI];(Reco - True) Muom Mom. (GeV) [SSRI]",  20, 0., 10., 5000, -10., 10.)
        #hMomRes[g] = ROOT.TH2D( "momRES_%s" % geom, ";True Momentum (SSRI);Momentum Resolution (SSRI)", 40, 0., 20., 100, 0., 1. )
        hMomRes[g] = ROOT.TH2D( "momRES_%s" % geom, ";#mu KE (GeV) [SSRI];Fractional Residual", 20, 0., 5., 100, -1.,1. )

        hVXY[g] = ROOT.TH2D( "denomvxy_%s" % geom, ";Vertex X;vertex Y", 40, -350., 350., 40, -400., 400. )
        hVYZ[g] = ROOT.TH2D( "denomvyz_%s" % geom, ";Vertex Z;vertex Y", 100, 0., 2000., 40, -400., 400. )
        hMuonBirthYZ[g] = ROOT.TH2D( "birthyz_%s" % geom, ";#mu Entry Z (cm); #mu Entry Y (cm)", 100, 500., 1500., 30, -400., 200. )
        hMuonBirthXZ[g] = ROOT.TH2D( "birthxz_%s" % geom, ";#mu Entry Z (cm);#mu Entry X (cm)", 100, 500., 1500., 40, -400., 400. )
        hMuonDeathYZ[g] = ROOT.TH2D( "deathyz_%s" % geom, ";#mu Endpoint Z (cm); #mu Endpoint Y (cm)", 100, 500., 1500., 30, -400., 200. )
        hMuonDeathXZ[g] = ROOT.TH2D( "deathxz_%s" % geom, ";#mu Endpoint Z (cm);#mu Endpoint X (cm)", 100, 500., 1500., 1600, -400., 400. )
        hDenom[g] = ROOT.TH2D( "denom_%s" % geom, ";Muon kinetic energy (GeV);Muon angle (degrees)", 30, 0., 6., 18, 0., 180. )
        hDenomQ2[g] = ROOT.TH2D( "denomQ2_%s" % geom, ";Neutrino energy (GeV);Q^{2} (GeV^{2})", 25, 0., 5., 25, 0., 5. )

        hDenomX[g] = ROOT.TH2D( "denomX_%s" % geom, ";Muon kinetic energy (GeV);Vertex X (cm)", 30, 0., 6., 60, -300., 300. )
        hNumX[g] = ROOT.TH2D( "numX_%s" % geom, ";Muon kinetic energy (GeV);Vertex X (cm)", 30, 0., 6., 60, -300., 300. )
        for r,rec in enumerate(recs):
            hNum[g][r] = ROOT.TH2D( "num_%s_%s" % (geom, rec), ";Muon kinetic energy (GeV);Muon angle (degrees)", 30, 0., 6., 18, 0., 180. )
            hNumLen[g][r] = ROOT.TH2D( "numlen_%s_%s" % (geom, rec), ";Muon kinetic energy (GeV);Muon track length", 30, 0., 6., 20, 0., 2000. )
            hNumQ2[g][r] = ROOT.TH2D( "numQ2_%s_%s" % (geom, rec), ";Neutrino energy (GeV);Q^{2} (GeV^{2})", 25, 0., 5., 25, 0., 5. )
        for t, ty in enumerate(nutype):
            hSignDist[g][t] = ROOT.TH1D( "signdist_%s_%s" % (geom, ty), ";Signed Distance (cm)", 100, -300., 300. )
            #hSignDistMom[g][t] = ROOT.TH2D( "signdistmom_%s_%s" % (geom, ty), ";True Muon momentum (GeV);Signed Distance (cm)", 20, 0., 20., 100, -300., 300. )
            hSignDistMom[g][t] = ROOT.TH2D( "signdistmom_%s_%s" % (geom, ty), ";Signed Distance (cm);Muon momentum (GeV) [LAr-exit]", 100, -300., 300., 60, 0., 6. )
        print("Looping for %s" % geom)
        loop( tree, hNum[g], hNumLen[g], hDenom[g], hNumX[g], hDenomX[g], hNumQ2[g], hDenomQ2[g], hVXY[g], hVYZ[g], hMuonBirthXZ[g], hMuonBirthYZ[g], hMuonDeathXZ[g], hMuonDeathYZ[g], hAngle[g], hED[g], hED2[g], hEDTrk[g], hScint[g], hMomAC[g], hMomSSRI[g], hMomCON[g], hMomRes[g], hMomRes2[g], hLArFactor[g], hSignDist[g], hSignDistMom[g] )

        #out = ROOT.TFile( "%s/filled_fifty.root" % topdir, "RECREATE" )
    out = ROOT.TFile( "%s/dumpSSRI-v5-histos.root" % topdir, "RECREATE" )
    for g,geom in enumerate(geometries):
        hDenom[g].Write()
        hLArFactor[g].Write()
        hDenomQ2[g].Write()
        hNumX[g].Write()
        hDenomX[g].Write()
        hVXY[g].Write()
        hVYZ[g].Write()
        hAngle[g].Write()
        hED[g].Write()
        hED2[g].Write()
        hEDTrk[g].Write()
        hScint[g].Write()
        hMomAC[g].Write()
        hMomSSRI[g].Write()
        hMomCON[g].Write()
        hMomRes[g].Write()
        hMomRes2[g].Write()
        hMuonBirthXZ[g].Write()
        hMuonBirthYZ[g].Write()
        hMuonDeathXZ[g].Write()
        hMuonDeathYZ[g].Write()

        for r,rec in enumerate(recs):
            hNum[g][r].Write()
            hNumLen[g][r].Write()
            hNumQ2[g][r].Write()


        for t, ty in enumerate(nutype):
            hSignDist[g][t].Write()
            hSignDistMom[g][t].Write()

