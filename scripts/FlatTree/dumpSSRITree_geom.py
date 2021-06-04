#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array
import glob

muon_mass = 105.6583755
nmax = -1
maxfiles = 9999

def ResetVariables():

  t_Reac = "EMPTY"

  ## initialize output variables
  t_vtx = [-999]*3
  t_p3lep = [-999]*3
  t_GlobalDeath = [-999]*3

  # In TMS
  t_TMSDeath = [-999]*3
  t_TMSBirth = [-999]*3

  # Exiting the LAr
  t_LAr_ExitPt = [-999]*3
  t_LAr_ExitMom = [-999]*3
  t_TMSLep = [-999]*3

  t_lepPdg[0] = -999
  t_lepKE[0] = -999
  t_lepE[0] = -999
  t_TMSKE[0] = -999
  t_Ev[0] = -999
  t_PDGv[0] = -999

  t_muonExitKE[0] = -999
  t_muonReco[0] = -999
  t_muonStart[0] = -999

  t_nFS = 0
  t_fsPdg = [-999]*100
  t_fsPx = [-999]*100
  t_fsPy = [-999]*100
  t_fsPz = [-999]*100
  t_fsE = [-999]*100

  t_muScintLen[0] = 0.0;
  t_muLArLen[0] = 0.0;
  t_muScintEnergy[0] = 0.0;

  xpt.clear()
  ypt.clear()
  zpt.clear()
  deposit.clear()


def loop( events, dspt, tgeo, tout ):

    offset = [ 0., 0., 0. ]

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    print "Starting loop over %d entries" % N
    ient = 0
    for evt in dspt:

        # Check if we're above maximum requested
        if ient > nmax: 
          break

        # Get the events
        events.GetEntry(ient)

        # See edep-sim/src/kinem/EDepSimRooTrakcerKinematicsGenerator
        t_Ev[0] = evt.StdHepP4[0][3]
        t_PDGv[0] = evt.StdHepPdg[3]

        if ient % 1000 == 0:
            print "Event %d of %d..." % (ient,N)


        # Loop over the primary vertices in the event
        for ivtx,vertex in enumerate(event.Primaries):

            # Reset the variables
            ResetVariables()

            # Save the event number
            t_ievt[0] = ient

            # Get the reaction code
            t_Reac = vertex.GetReaction()

            # Set the vertex location
            for i in range(3): t_vtx[i] = vertex.GetPosition()[i] / 10. - offset[i] # cm

            # Reset lepton trajectory
            ileptraj = -1
            # Number of final state particle
            nfsp = 0

            # Loop over the particles at the vertex
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.GetMomentum()[3]
                # Get the momentum (px^2+py^2+pz^2)
                p = (particle.GetMomentum()[0]**2 + particle.GetMomentum()[1]**2 + particle.GetMomentum()[2]**2)**0.5
                m = (e**2 - p**2)**0.5
                pdg = particle.GetPDGCode()
                t_fsPdg[nfsp] = pdg
                t_fsPx[nfsp] = particle.GetMomentum()[0]
                t_fsPy[nfsp] = particle.GetMomentum()[1]
                t_fsPz[nfsp] = particle.GetMomentum()[2]
                t_fsE[nfsp] = e
                nfsp += 1

                # Look for the lepton (dicey, need match?)
                if abs(pdg) in [11,12,13,14]:
                    ileptraj = particle.GetTrackId()
                    t_lepPdg[0] = pdg
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.GetMomentum()[i]
                    t_lepKE[0] = e - m
                    t_lepE[0] = e

            # Finish the loop
            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp

            # Only interested in muons: move to next vertex if not muon
            if abs(t_lepPdg[0]) != 13: continue

            # Check if this trajectory exits the LAr, and is in the TMS
            inTMS = False

            # The lepton trajectory
            leptraj = event.Trajectories[ileptraj]

            # Can loop over each trajectory
            #for traj in event.Trajectories:

            # Now find which volume the first hit of the lepton trajectory is
            startpt = leptraj.Points[0].GetPosition()
            startnode = tgeo.FindNode( startpt.X(), startpt.Y(), startpt.Z() )
            vertexname = startnode.GetName()
            if ( "volTPCActive_" in vertexname or 
                 "volPixelPlane_" in vertexname or 
                 "volLAr_" in vertexname or
                 "volArCLight_PV" in vertexname ):
              t_muonStart[0] = 1
            elif ("TMS" in vertexname or
                  "modulelayer" in vertexname):
              t_muonStart[0] = 2
            else:
              t_muonStart[0] = 0

            # Now loop over the points in the trajectory
            for p in leptraj.Points:
                pt = p.GetPosition()
                node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                volName = node.GetName()

                pPos = ROOT.TVector3( pt.X()/10. - offset[0], pt.Y()/10. - offset[1], pt.Z()/10. - offset[2] )
                # Check where particle is, update if it's in active LAr
                if ( "volTPCActive_" in volName or 
                     "volPixelPlane_" in volName or 
                     "volLAr_" in volName or
                     "volArCLight_PV" in volName ):
                    t_LAr_ExitPt[0] = pt.X() / 10. - offset[0]
                    t_LAr_ExitPt[1] = pt.Y() / 10. - offset[1]
                    t_LAr_ExitPt[2] = pt.Z() / 10. - offset[2]
                    t_LAr_ExitMom[0] = p.GetMomentum().x()
                    t_LAr_ExitMom[1] = p.GetMomentum().y()
                    t_LAr_ExitMom[2] = p.GetMomentum().z()
                    t_muonExitKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass

                # Else the particle has exited
                else:
                    # Check if it is in the TMS
                    if ("TMS_" in volName or 
                        "modulelayer" in volName) and not inTMS:
                        t_TMSKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                        inTMS = True

            # Get end point
            endpt = leptraj.Points[-1].GetPosition()

            # Find the node correpsonding to the end point
            node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )
            # Save the death point
            t_GlobalDeath[0] = endpt.X()/10. - offset[0]
            t_GlobalDeath[1] = endpt.Y()/10. - offset[1]
            t_GlobalDeath[2] = endpt.Z()/10. - offset[2]

            # Check which volume the muon dies in and save containment
            endVolName = node.GetName()
            if "TPCActive" in endVolName: t_muonReco[0] = 1 # contained
            # Updated name in new geom
            elif "TMS" in endVolName: t_muonReco[0] = 2 # Scintillator stopper
            else: t_muonReco[0] = 0 # endpoint not in active material

            # look for muon hits in the ArgonCube
            arhits = []
            for key in event.SegmentDetectors:
                # Updated name in new geom
                if key.first == "TPCActive_shape":
                    arhits += key.second
                    
            ar_muon_hits = []
            for idx, hit in enumerate(arhits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                #if traj.GetParentId() == -1 and abs(traj.GetPDGCode()) == 13:
                ar_muon_hits.append(hit)

            for idx, hit in enumerate(ar_muon_hits):
                hStart = ROOT.TVector3( hit.GetStart()[0]/10.-offset[0], hit.GetStart()[1]/10.-offset[1], hit.GetStart()[2]/10.-offset[2] )
                xpt.push_back(hStart.x())
                ypt.push_back(hStart.y())
                zpt.push_back(hStart.z())
                
            # Use vertex and LAr endpoint for first pass
            # Maybe not the best idea...
            hArStart=ROOT.TVector3(ar_muon_hits[0].GetStart()[0]/10.-offset[0], ar_muon_hits[0].GetStart()[1]/10.-offset[1], ar_muon_hits[0].GetStart()[2]/10.-offset[2])
            hArEnd=ROOT.TVector3(ar_muon_hits[-1].GetStart()[0]/10.-offset[0], ar_muon_hits[-1].GetStart()[1]/10.-offset[1], ar_muon_hits[-1].GetStart()[2]/10.-offset[2])

            ar_dh = (hArEnd-hArStart).Mag()
            ar_dz = hArEnd.z()-hArStart.z()
                
            ar_trk_length_gcm2 = 1.4 * ar_dh

            # Muon hits in the scintillator volumes
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "volTMS":
                    hits += key.second

            muon_hits = []
            for idx, hit in enumerate(hits):
                tid = hit.Contrib[0]
                traj = event.Trajectories[tid]
                # Check they are fundamental
                #if traj.GetParentId() == -1 and abs(traj.GetPDGCode()) == 13:
                muon_hits.append(hit)

            if muon_hits:
              # If the track enters the TMS, it has passed through the passive material and one layer of 1.5cm steel
              # Accounting for the angle, add up the track length
              hMuonStart = muon_hits[0].GetStart()
              t_TMSBirth[0] = hMuonStart[0]/10.-offset[0]
              t_TMSBirth[1] = hMuonStart[1]/10.-offset[1]
              t_TMSBirth[2] = hMuonStart[2]/10.-offset[2]

              # Track length
              trk_length_gcm2 = 0.
              # Total deposited energy
              de = 0.

              # Reconstruct the exit point in the LAr
              LArExit = ROOT.TVector3(t_LAr_ExitPt[0], t_LAr_ExitPt[1], t_LAr_ExitPt[2])
              # Add in the extra length from traversing the dead region
              costh = (hMuonStart.Z()-LArExit.Z())/(hMuonStart-LArExit).Mag()

            hPrev = None
            hFinal = None
            hPrevTime = 0
            total_planes = 0
            prevlayerno = -1

            for idx, hit in enumerate(muon_hits):
                de += hit.GetEnergyDeposit()

                hStart = ROOT.TVector3( hit.GetStart()[0]/10.-offset[0], hit.GetStart()[1]/10.-offset[1], hit.GetStart()[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.GetStop()[0]/10.-offset[0], hit.GetStop()[1]/10.-offset[1], hit.GetStop()[2]/10.-offset[2] )
                Time = hit.GetStart()[3]

                xpt.push_back(hStart.x())
                ypt.push_back(hStart.y())
                zpt.push_back(hStart.z())
                deposit.push_back(hit.GetEnergyDeposit())

                # Check if this point is in the same bar as the one before it
                # Move on if it is
                same = tgeo.IsSameLocation( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2])

                tempnode = tgeo.FindNode( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2])
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

                # In x, gaps at -177 to -175 roughly
                #               +175 to +177
                # In y, gaps are NOWHERE
                # In z, gaps are

                if hPrev is not None:
                    gap = hStart.z() - hPrev.z()
                    # correct for cosine of the incident angle
                    dh = (hStart-hPrev).Mag()
                    dz = hStart.z()-hPrev.z()


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

            t_TMSDeath[0] = hFinal.x()
            t_TMSDeath[1] = hFinal.y()
            t_TMSDeath[2] = hFinal.z()

            # Take the directionality into account when adding the track length
            t_muScintLen[0] = trk_length_gcm2
            t_muLArLen[0] = ar_trk_length_gcm2 
            t_muScintEnergy[0] = de
            
            tout.Fill()

        ient += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--nmax', help='Maximum number of events', default="99999999999999")
    parser.add_option('--maxfiles', help='Maximum number of events', default="99999")
    #parser.add_option('--first_run', type=int, help='First run number', default=0)
    #parser.add_option('--last_run', type=int, help='Last run number', default=0)
    #parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    #parser.add_option('--grid', action='store_true', help='grid mode', default=False)

    (args, dummy) = parser.parse_args()

    nmax = int(args.nmax)
    maxfiles = int(args.maxfiles)

    # Create the output tree
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )

    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_Ev = array('f', [0.])
    tout.Branch('Ev',t_Ev,'Ev/F')
    t_PDGv = array('i', [0])
    tout.Branch('PDGv',t_PDGv,'PDGv/I')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')

    t_GlobalDeath = array('f',3*[0.0])
    tout.Branch('GlobalDeath',t_GlobalDeath,'GlobalDeath[3]/F')
    t_TMSDeath = array('f',3*[0.0])
    tout.Branch('TMSDeath',t_TMSDeath,'TMSDeath[3]/F')
    t_TMSBirth = array('f',3*[0.0])
    tout.Branch('TMSBirth',t_TMSBirth,'TMSBirth[3]/F')
    t_TMSLep = array('f',3*[0.0])
    tout.Branch('TMSLep',t_TMSLep,'leTMS[3]/F')

    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_TMSKE = array('f',[0])
    tout.Branch('TMSKE',t_TMSKE,'TMSKE/F')
    t_lepE = array('f',[0])
    tout.Branch('lepE',t_lepE,'lepE/F')

    t_LAr_ExitPt = array('f',3*[0.0])
    tout.Branch('LAr_ExitPt',t_LAr_ExitPt,'LAr_ExitPt[3]/F')
    t_LAr_ExitMom = array('f',3*[0.0])
    tout.Branch('LAr_ExitMom',t_LAr_ExitMom,'LAr_ExitMom[3]/F')

    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muonStart = array('i',[0])
    tout.Branch('muonStart',t_muonStart,'muonStart/I')

    # Track-length quantities
    t_TMS_TrackLength = array('f',[0])
    tout.Branch('TMS_TrackLength', t_TMS_TrackLength, 'TMS_TrackLength/F')
    t_LAr_TrackLength = array('f',[0])
    tout.Branch('LAr_TrackLength',t_LAr_TrackLength,'LAr_TrackLength/F')

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

    t_Reac = "";
    tout.Branch('reac', t_Reac, 'reac/C')

    xpt = ROOT.std.vector('float')()
    ypt = ROOT.std.vector('float')()
    zpt = ROOT.std.vector('float')()
    deposit = ROOT.std.vector('float')()
    tout.Branch('xpt', xpt)
    tout.Branch('ypt', ypt)
    tout.Branch('zpt', zpt)
    tout.Branch('deposit', deposit)
    
    loaded = False
    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )
    dspt = ROOT.TChain( "DetSimPassThru/gRooTracker", "other thing" )

    #neutrino = "neutrino"
    #if args.rhc:
        #neutrino = "antineutrino"

    #print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    #for run in range( args.first_run, args.last_run+1 ):
        #if args.grid:
            #fname = "%s/edep.%d.root" % (args.topdir,run)
        #else:
            #fname = "%s/%s.%d.edepsim.root" % (args.topdir, neutrino, run)
        #print fname
#
        #fname = fname.replace("/pnfs","root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr")
        # see if it is an OK file
        #if not os.access( fname, os.R_OK ):
            #print "Can't access file: %s" % fname
            #continue
        #tf = ROOT.TFile.Open( fname )
        #if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            #print "File is crap: %s" % fname
            #continue
#
        #if not loaded:
            #loaded = True
            #tf.MakeProject("EDepSimEvents","*","RECREATE++")
##
        # add it to the tchain
        #events.Add( fname )
        #dspt.Add( fname )

        #if tgeo is None: # first OK file, get geometry
            #tgeo = tf.Get("EDepSimGeometry")
        #tf.Close() # done with this one

    nfiles = 0
    for run in glob.glob(args.topdir+"/edep*.root"):
      print nfiles
      print maxfiles
      if nfiles > maxfiles: 
        break

      fname = run
      print "Adding "+fname+" to TChain..."

      #fname = fname.replace("/pnfs","root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr")
      # see if it is an OK file
      if not os.access( fname, os.R_OK ):
        print "Can't access file: %s" % fname
        continue

      tf = ROOT.TFile.Open( fname )
      if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
        print "File is crap: %s" % fname
        continue

      if not loaded:
        loaded = True
        tf.MakeProject("EDepSimEvents","*","RECREATE++")

      # add it to the tchain
      events.Add( fname )
      dspt.Add( fname )

      # first OK file, get geometry
      if tgeo is None:
        tgeo = tf.Get("EDepSimGeometry")
      tf.Close() # done with this one

      nfiles += 1

    print "Running on ",nfiles," root files..."
    loop( events, dspt, tgeo, tout )

    fout.cd()
    tout.Write()

