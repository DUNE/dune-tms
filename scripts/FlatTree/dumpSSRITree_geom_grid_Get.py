#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
import glob

muon_mass = 105.6583755
extra_trk_length = 24.35 # g/cm^2 from Mike K.
nmax = -1
maxfiles = 9999

def loop( events, dspt, tgeo, tout ):

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    print "Starting loop over %d entries" % N
    ient = 0
    for evt in dspt:

        t_Ev[0] = evt.StdHepP4[3] # neutrino is first 4-vector, so 3 is always Ev

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)

        if ient > nmax: 
          break

        events.GetEntry(ient)

        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            t_ievt[0] = ient
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            t_lepBirth[0]=0.0; t_lepBirth[1]=0.0; t_lepBirth[2]=0.0;
            # In SSRI
            t_muonDeath[0]=0.0; t_muonDeath[1]=0.0; t_muonDeath[2]=0.0;
            t_muonBirth[0]=0.0; t_muonBirth[1]=0.0; t_muonBirth[2]=0.0;

            t_lepRMMS[0]=0.0; t_lepRMMS[1]=0.0; t_lepRMMS[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_rmmsKE[0] = -1.
            t_lepE[0] = 0.
            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonExitKE[0] = 0.0
            t_muonReco[0] = -1;
            t_muScintLen[0] = 0.0;
            t_muLArLen[0] = 0.0;
            t_muScintEnergy[0] = 0.0;

            xpt.clear()
            ypt.clear()
            zpt.clear()
            deposit.clear()

            #pxpt.clear()
            #pypt.clear()
            #pzpt.clear()
            #px.clear()
            #py.clear()
            #pz.clear()

            # now ID numucc
            #t_Reac=vertex.Reaction

            # set the vertex location for output
            for i in range(3): t_vtx[i] = vertex.GetPosition()[i]

            ileptraj = -1
            nfsp = 0

            # get the lepton kinematics from the edepsim file
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.GetMomentum()[3]
                # Get the momentum (px^2+py^2+pz^2)
                p = (particle.GetMomentum()[0]**2 + particle.GetMomentum()[1]**2 + particle.GetMomentum()[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                m = (e**2 - p**2)**0.5
                t_fsPdg[nfsp] = particle.GetPDGCode()
                t_fsPx[nfsp] = particle.GetMomentum()[0]
                t_fsPy[nfsp] = particle.GetMomentum()[1]
                t_fsPz[nfsp] = particle.GetMomentum()[2]
                t_fsE[nfsp] = e
                nfsp += 1
                pdg = particle.GetPDGCode()

                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.GetTrackId()
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.GetMomentum()[i]
                    t_lepKE[0] = e - m
                    t_lepE[0] = e

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp


            if abs(t_lepPdg[0]) != 13: continue

            # If there is a muon, determine how to reconstruct its momentum and charge
            exit = False
            inrmms = False

            leptraj = event.Trajectories[ileptraj]

            pPrev = None
            pointval = 0
            for p in leptraj.Points:
                pt = p.GetPosition()
                node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                volName = node.GetName()
                active = False

                pPos = ROOT.TVector3( pt.X(), pt.Y(), pt.Z() )
                #print "point ",pointval,pPos.z(),volName

                # The z-information is contained in thinlayerposition0, thinlayervol

                if "LAr" in volName or "PixelPlane" in volName or "sPlane" in volName: # in active volume, update exit points
                    t_muonExitPt[0] = pt.X()
                    t_muonExitPt[1] = pt.Y()
                    t_muonExitPt[2] = pt.Z()
                    t_muonExitMom[0] = p.GetMomentum().x()
                    t_muonExitMom[1] = p.GetMomentum().y()
                    t_muonExitMom[2] = p.GetMomentum().z()
                    
                    t_muonExitKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                else:
                    if not exit:
                        t_muonExitPt[0] = pt.X()
                        t_muonExitPt[1] = pt.Y()
                        t_muonExitPt[2] = pt.Z()
                        exit = True

                    # Check if it is in the RMMS
                    if ("RMMS" in volName or "modulelayer" in volName) and not inrmms:
                        t_rmmsKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                        inrmms = True
                        print pPos.Z()
                pPrev = pPos
                pointval += 1
                # Save up the momentum
                #px.push_back(p.GetMomentum().x())
                #py.push_back(p.GetMomentum().y())
                #pz.push_back(p.GetMomentum().z())
                # And the position
                #pxpt.push_back(pPos.x())
                #pypt.push_back(pPos.y())
                #pzpt.push_back(pPos.z())

            endpt = leptraj.Points[-1].GetPosition()
            startpt = leptraj.Points[0].GetPosition()

            node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

            t_lepBirth[0] = startpt.X()
            t_lepBirth[1] = startpt.Y()
            t_lepBirth[2] = startpt.Z()

            t_lepDeath[0] = endpt.X()
            t_lepDeath[1] = endpt.Y()
            t_lepDeath[2] = endpt.Z()

            endVolName = node.GetName()
            if "LArActive" in endVolName: t_muonReco[0] = 1 # contained
            elif "RMMS" in endVolName: t_muonReco[0] = 2 # Scintillator stopper
            else: t_muonReco[0] = 0 # endpoint not in active material

            # look for muon hits in the ArgonCube
            arhits = []
            for key in event.SegmentDetectors:
                if key.first == "ArgonCube":
                    arhits += key.second
                    
            ar_muon_hits = []
            for idx, hit in enumerate(arhits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                if traj.ParentId == -1 and abs(traj.GetPDGCode()) == 13:
                    ar_muon_hits.append(hit)

            #if len(ar_muon_hits) < 3: continue

            ar_trk_length_gcm2 = 0.
            if len(ar_muon_hits) > 2: 
              for idx, hit in enumerate(ar_muon_hits):
                  hStart = ROOT.TVector3( hit.Start[0], hit.Start[1], hit.Start[2] )
                  xpt.push_back(hStart.x())
                  ypt.push_back(hStart.y())
                  zpt.push_back(hStart.z())
                
            # just use vertex and LAr endpoint for first pass
              hArStart=ROOT.TVector3(ar_muon_hits[0].Start[0], ar_muon_hits[0].Start[1], ar_muon_hits[0].Start[2])
              hArEnd=ROOT.TVector3(ar_muon_hits[-1].Start[0], ar_muon_hits[-1].Start[1], ar_muon_hits[-1].Start[2])

              ar_dh = (hArEnd-hArStart).Mag()
              ar_dz = hArEnd.z()-hArStart.z()
                
              ar_trk_length_gcm2 = 1.4 * ar_dh

            #-------------------------------------------------------
            # look for muon hits in the scintillator
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "rmmsvol":
                    hits += key.second

            muon_hits = []
            for idx, hit in enumerate(hits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                if traj.ParentId == -1 and abs(traj.GetPDGCode()) == 13:
                    muon_hits.append(hit)

            # Need at least 2 hits in the TMS to get start and end
            #if len(muon_hits) < 2: continue

            # Track length
            trk_length_gcm2 = 0.
            # Deposited energy
            de = 0.

            # If the track enters the TMS, it has passed through the passive material and one layer of 1.5cm steel
            # Accounting for the angle, add up the track length
            if len(muon_hits) > 0:
              hMuonStart = muon_hits[0].Start.Vect()
              t_muonBirth[0] = hMuonStart[0]
              t_muonBirth[1] = hMuonStart[1]
              t_muonBirth[2] = hMuonStart[2]
              # Reconstruct the exit point in the LAr
              LArExit = ROOT.TVector3(t_muonExitPt[0], t_muonExitPt[1], t_muonExitPt[2])
              # Add in the extra length from traversing the dead region
              costh = (hMuonStart.Z()-LArExit.Z())/(hMuonStart-LArExit).Mag()
              trk_length_gcm2 += extra_trk_length/costh
              # And the extra length from the first thin layer of steel and scintillator
              # Should the scinillator really be added?
              #trk_length_gcm2 += (1.05+1.5*7.85)/costh
              trk_length_gcm2 += (1.5*7.85)/costh

            hPrev = None
            hFinal = None
            hPrevTime = 0
            total_planes = 0
            prevlayerno = -1

            for idx, hit in enumerate(muon_hits):
                de += hit.EnergyDeposit

                hStart = ROOT.TVector3( hit.Start[0], hit.Start[1], hit.Start[2])
                hStop = ROOT.TVector3( hit.Stop[0], hit.Stop[1], hit.Stop[2])
                Time = hit.Start[3]

                xpt.push_back(hStart.x())
                ypt.push_back(hStart.y())
                zpt.push_back(hStart.z())
                deposit.push_back(hit.EnergyDeposit)

                # Check if this point is in the same bar as the one before it
                # Move on if it is
                same = tgeo.IsSameLocation( hit.Start[0], hit.Start[1], hit.Start[2])

                tempnode = tgeo.FindNode( hit.Start[0], hit.Start[1], hit.Start[2])
                nav = tgeo.GetCurrentNavigator()
                topnode_name = nav.GetCurrentNode().GetName()

                # Z enumerated volume (each bar)
                while "modulelayervol_PV" not in topnode_name and "volWorld" not in topnode_name:
                  nav.CdUp()
                  topnode_name = nav.GetCurrentNode().GetName()

                # No SSRI hits (cd through the geom and found nothing)
                if "volWorld" in topnode_name:
                  continue

                layerno = nav.GetCurrentNode().GetNumber()

                nav.GetCurrentNode()

                # Ignore double hits in a bar for track length calculation
                if layerno == prevlayerno or same == True:
                  continue

                # Don't want delayed hits (only look forwards)
                if layerno < prevlayerno:
                  continue

                Consecutive = True
                nlayers = layerno-prevlayerno
                if nlayers != 1:
                  Consecutive = False
                  #print "not consecutive"
                  #if hPrev is not None:
                    #hPrev.Print()
                  #hStart.Print()

                if hPrev is not None:
                    # thin layer is 3cm air + 1cm scint + 1.5cm steel = 5.5cm pitch
                    # -> means 3cm+1.5cm = 4.5cm between hits
                    # Runs from 733 to 949

                    # thick layer is 3cm air + 1cm scint + 4cm steel = 8cm pitch
                    # -> means 3cm+4cm = 7cm between hits

                    # There is also the transition layer
                    # 1.5cm steel + 4cm steel stuck together, with air on either side (each 1.5cm) 
                    # -> 1.5+4cm+1.5+1.5 = 8.5cm between hits
                    # find 8.5cm between hits

                    # if "gap" isn't a multiple of that, then something has gone off the rails and we want to stop
                    gap = hStart.z() - hPrev.z()
                    # correct for cosine of the incident angle
                    dh = (hStart-hPrev).Mag()
                    dz = hStart.z()-hPrev.z()

                    #print "prevlayerno",prevlayerno
                    #print "layerno",layerno
                    #print "gap",gap

                    #Found change at 730.8-736.3
                    #Previous distance: 730.8
                    #New distance: 5.5

                    #Found change at 939.8-949.3
                    #Previous distance: 5.5
                    #New distance: 9.5

                    #Found change at 949.3-957.3
                    #Previous distance: 9.5
                    #New distance: 8

                    # Minimum gap between hits is end of scintillator bar n to beginning of scintillator bar n+1: i.e. 1.5cm air + 1.5cm steel + 1.5cm air -> 4.5cm
                    if gap < 4.5:
                      prevlayerno = layerno
                      hPrev = hStart
                      hPrevTime = Time
                      continue

                    # Make a cut on the time difference to avoid delayed hits
                    # Only when the xhit is not in the dead region
                    # Time between hits can be large if the particle propagates in a region without scintillator bars
                    timediff = Time-hPrevTime
                    if timediff > 2 and nlayers < 2:
                      prevlayerno = layerno
                      hPrev = hStart
                      hPrevTime = Time
                      continue

                    # Increment the total planes traversed
                    total_planes += nlayers

                    # Account for possibility that we've had no hits whilst transitioning over more than one plane (e.g. traversed region where there is iron but no scint)
                    # Requires a hit in the thin or transition layer
                    if nlayers > 1:

                      # First make sure the trajectory didn't exit the physical boundary of the detector
                      # The gap is around 170 to 180 cm in x
                      if abs(hStart.x()) > 200:
                        # If it exited, it definitely did pass through a layer of scintillator bar to make this hit
                        trk_length_gcm2 += 1.05

                        hPrev = hStart
                        prevlayerno = layerno
                        hPrevTime = Time
                        #print "trklen",trk_length_gcm2
                        continue

                      # The current layer is in thick or transition region
                      # The previous layer was in the thin region
                      if layerno >= 39 and prevlayerno < 39:
                        # Count up layers traversed in thin region (up to layer 0,...,38
                        nthins = 37-prevlayerno
                        # Check if we've traversed the transition layer (between 38 and 39)
                        if layerno > 38:
                          ntrans = 1
                        else:
                          ntrans = 0

                        nthick = layerno-39
                        #print "thin,trans,thick=",nthins, ntrans, nthick

                        # Since we had no hits in the scintiallator bars here, we only went through iron
                        # In thin region, that's 1.5cm iron
                        trk_length_gcm2 += (1.5*7.85) * nthins * dh/dz
                        # Account for transition region, where it's a 5.5cm iron
                        trk_length_gcm2 += ((1.5+4.0)*7.85) * ntrans * dh/dz
                        # Account for thick layers, where it's 4.5cm iron
                        trk_length_gcm2 += (4.0*7.85) * nthick * dh/dz

                        #print "trklen",trk_length_gcm2

                      # If hits in adjacent bars, it went through scintillator and iron
                      # If not, it did not pass through scintillator, only iron
                      else:

                        if layerno > 39 and prevlayerno > 39:
                          nthick = layerno-prevlayerno
                          trk_length_gcm2 += (4.0*7.85) * nthick * dh/dz
                        elif layerno < 39 and prevlayerno < 39:
                          nthin = layerno-prevlayerno
                          trk_length_gcm2 += (1.5*7.85) * nthin * dh/dz

                        #print "trklen",trk_length_gcm2


                    # else we just calculate is as per usual
                    else:

                    # Anything before 39 is thin module
                      if layerno < 39:
                      #print "thin ",gap
                      # Get the gap right from the geometry
                        nplanes = int(gap/5.5)
                      # Should be equal to layerno-prevlayerno

                        trk_length_gcm2 += (1.05 + 1.5*7.85) * nlayers * dh/dz

                    # Layer number 39 is transition layer
                      elif layerno == 39:
                      #print "trans ", gap
                        nplanes = int(gap/9.5)

                      # A hit here definitely means it passed through scintillator and transition
                        trk_length_gcm2 += (1.05+(1.5+4.0)*7.85) * nlayers * dh/dz

                    # Anything after 39 is thick module
                      else:
                      #print "thick ", gap
                        nplanes = int(gap/8.0)
                        trk_length_gcm2 += (1.05+4.0*7.85) * nlayers * dh/dz

                      #print "trklen",trk_length_gcm2

                      '''
                      if abs(nplanes - nlayers)>1:
                        print "***"
                        #print "EHH", nplanes, gap, timediff
                        print "EHH", nplanes, gap
                        print nlayers, layerno, prevlayerno
                        hPrev.Print()
                        hStart.Print()
                        sys.exit(-1)
                      '''

                      '''
                        nplanes = int(gap/5.5)
                        # if nplanes isn't an integer, then we've reached the end of the track
                        # and there is some other hit that is super far away typically that is messing up the length

                        total_planes += nplanes
                        # Sometimes the track goes out of scintillator, still propagates through the iron and air, then comes back into the scintillator
                        # This leads to a larger gap, but it only goes through the iron, so modify track length
                        # The dh/dz modification (costheta) will also be borked (have no info on bend), so assume straight line
                        if gap > 5.5 and nplanes > 1:
                          trk_length_gcm2 += (1.5*7.85) * nplanes * dh/dz
                        else:
                          # we must have gone through 1cm scint + 1.5cm steel + 3cm air, divide by cos(theta)
                          trk_length_gcm2 += (1.05 + 1.5*7.85) * nplanes * dh/dz

                    # The transition layer gap
                    # Between 939.8+1 to 949.3
                    #elif (hStart.z() - 949.3) < 0.1: # the transition module, the gap is 9.5

                    # Scintillator layer 39/40
                    elif (hStart.z() < 949.3+1.5): # Scintillator strip goes from 949.3 to 949.3+1 (add in 1.5 for good measure)
                        print layerno
                        if abs(gap-9.5) < 1.E-3: # Gap is 9.5cm
                            trk_length_gcm2 += (1.05 + (4.+1.5)*7.85) * dh/dz # Gone through 1 cm scintillator, 4cm steel, 1.5cm steel, 3cm air
                            total_planes += 1 # Add in the plane

                    else: #"thick layer region, gap should be 8cm
                        nplanes = int(gap/8.)
                        total_planes += nplanes

                        if gap > 8 and nplanes > 1:
                          trk_length_gcm2 += (4.*7.85) * nplanes * dh/dz
                        else:
                          trk_length_gcm2 += (1.05 + 4.*7.85) * nplanes * dh/dz
                  '''

                # Update the previous layer number
                prevlayerno = layerno
                # Update start position
                hPrev = hStart
                # Update the time
                hPrevTime = Time
                # Update the final point
                hFinal = hStop

            if hFinal is not None:
              t_muonDeath[0] = hFinal.x()
              t_muonDeath[1] = hFinal.y()
              t_muonDeath[2] = hFinal.z()

            #t_muScintLen[0] = trk_length_gcm2 + extra_trk_length
            # Take the directionality into account when adding the track length
            t_muScintLen[0] = trk_length_gcm2 + extra_trk_length
            t_muLArLen[0] = ar_trk_length_gcm2 
            t_muScintEnergy[0] = de
            
            tout.Fill()
        ient += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--nmax', help='Maximum number of events', default="99999999999999")
    parser.add_option('--infile', help='EDepSim input file')

    (args, dummy) = parser.parse_args()

    nmax = int(args.nmax)
    input_file = args.infile
    print input_file
    if not input_file:
      print "Need input file"
      exit(-1)

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_Ev = array('f', [0.])
    tout.Branch('Ev',t_Ev,'Ev/F')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')

    t_lepBirth = array('f',3*[0.0])
    tout.Branch('lepBirth',t_lepBirth,'lepBirth[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')

    t_muonDeath = array('f',3*[0.0])
    tout.Branch('muonDeath',t_muonDeath,'muonDeath[3]/F')
    t_muonBirth = array('f',3*[0.0])
    tout.Branch('muonBirth',t_muonBirth,'muonBirth[3]/F')
    t_lepRMMS = array('f',3*[0.0])
    tout.Branch('lepRMMS',t_lepRMMS,'leRMMS[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_rmmsKE = array('f',[0])
    tout.Branch('rmmsKE',t_rmmsKE,'rmmsKE/F')
    t_lepE = array('f',[0])
    tout.Branch('lepE',t_lepE,'lepE/F')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muScintLen = array('f',[0])
    tout.Branch('muScintLen',t_muScintLen,'muScintLen/F')
    t_muLArLen = array('f',[0])
    tout.Branch('muLArLen',t_muLArLen,'muLArLen/F')
    t_muonExitKE = array('f',[0])
    tout.Branch('muonExitKE',t_muonExitKE,'muonExitKE/F')
    t_muScintEnergy = array('f',[0])
    tout.Branch('muScintEnergy',t_muScintEnergy,'muScintEnergy/F')
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')

    xpt = ROOT.std.vector('float')()
    ypt = ROOT.std.vector('float')()
    zpt = ROOT.std.vector('float')()
    deposit = ROOT.std.vector('float')()
    tout.Branch('xpt', xpt)
    tout.Branch('ypt', ypt)
    tout.Branch('zpt', zpt)
    tout.Branch('deposit', deposit)

    #px = ROOT.std.vector('float')()
    #py = ROOT.std.vector('float')()
    #pz = ROOT.std.vector('float')()

    #pxpt = ROOT.std.vector('float')()
    #pypt = ROOT.std.vector('float')()
    #pzpt = ROOT.std.vector('float')()

    #tout.Branch('px',px)
    #tout.Branch('py',py)
    #tout.Branch('pz',pz)
    #tout.Branch('pxpt',pxpt)
    #tout.Branch('pypt',pypt)
    #tout.Branch('pzpt',pzpt)
    
    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    fname = input_file
    print "Adding "+fname+" to TChain..."
    tgeo = None

    tf = ROOT.TFile.Open( fname )
    tf.MakeProject("EDepSimEvents","*","RECREATE++")

    # add it to the tchain
    events.Add( fname )
    dspt.Add( fname )

    tgeo = tf.Get("EDepSimGeometry")
    tf.Close() # done with this one

    loop( events, dspt, tgeo, tout )

    fout.cd()
    tout.Write()

