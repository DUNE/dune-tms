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

  t_Reac="EMPTY"

  ## initialize output variables
  t_vtx = [-999]*3
  t_p3lep = [-999]*3
  t_GlobalDeath = [-999]*3

  # In TMS
  t_TMSDeath = [-999]*3
  t_TMSBirth = [-999]*3

  # Exiting the LAr
  t_muonExitPt = [-999]*3
  t_muonExitMom = [-999]*3
  t_lepTMS = [-999]*3

  t_lepPdg[0] = -999
  t_lepKE[0] = -999
  t_TMSKE[0] = -999
  t_lepE[0] = -999
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


        if ient % 100 == 0:
          print "Event %d of %d..." % (ient,N)

        if ient > nmax: 
          break

        events.GetEntry(ient)

        # Loop over the primary vertices in the event
        for ivtx,vertex in enumerate(event.Primaries):

            # Reset the variables
            ResetVariables()
            t_ievt[0] = ient
            t_Ev[0] = evt.StdHepP4[3] # neutrino is first 4-vector, so 3 is always Ev

            # Get the reaction code
            t_Reac = vertex.GetReaction()

            # Set the vertex location
            for i in range(3): t_vtx[i] = vertex.GetPosition()[i]

            # Reset lepton trajectory
            ileptraj = -1
            # Number of final state particle
            nfsp = 0

            # Loop over the particles at the vertex
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.GetMomentum()[3]
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

            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp

            #if abs(t_lepPdg[0]) != 13: continue

            # If there is a muon, determine how to reconstruct its momentum and charge
            exit = False
            inTMS = False

            leptraj = event.Trajectories[ileptraj]

            pointval = 0

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

            # Now loop over the particles
            for p in leptraj.Points:
                pt = p.GetPosition()
                node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                volName = node.GetName()

                pPos = ROOT.TVector3( pt.X(), pt.Y(), pt.Z() )
                if ( "volTPCActive_" in volName or 
                     "volPixelPlane_" in volName or 
                     "volLAr_" in volName or
                     "volArCLight_PV" in volName ):
                    t_muonExitPt[0] = pt.X()
                    t_muonExitPt[1] = pt.Y()
                    t_muonExitPt[2] = pt.Z()
                    t_muonExitMom[0] = p.GetMomentum().x()
                    t_muonExitMom[1] = p.GetMomentum().y()
                    t_muonExitMom[2] = p.GetMomentum().z()
                    
                    t_muonExitKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                else:
                    exit = True

                    # Check if it is in the TMS
                    if ("TMS_" in volName or 
                        "modulelayer" in volName) and not inTMS:
                        t_TMSKE[0] = (p.GetMomentum().Mag2() + muon_mass*2)**0.5 - muon_mass
                        inTMS = True

                pointval += 1

            endpt = leptraj.Points[-1].GetPosition()

            node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

            t_GlobalDeath[0] = endpt.X()
            t_GlobalDeath[1] = endpt.Y()
            t_GlobalDeath[2] = endpt.Z()

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
                hStart = ROOT.TVector3( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2] )
                xpt.push_back(hStart.x())
                ypt.push_back(hStart.y())
                zpt.push_back(hStart.z())
                deposit.push_back(hit.GetEnergyDeposit())
                
            #-------------------------------------------------------
            # look for muon hits in the scintillator
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
              hMuonStart = muon_hits[0].GetStart()
              t_TMSBirth[0] = hMuonStart[0]
              t_TMSBirth[1] = hMuonStart[1]
              t_TMSBirth[2] = hMuonStart[2]

              for idx, hit in enumerate(muon_hits):

                hStart = ROOT.TVector3( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2] )
                hStop = ROOT.TVector3( hit.GetStop()[0], hit.GetStop()[1], hit.GetStop()[2] )
                Time = hit.GetStart()[3]

                xpt.push_back(hStart.x())
                ypt.push_back(hStart.y())
                zpt.push_back(hStart.z())
                deposit.push_back(hit.GetEnergyDeposit())

              hit = muon_hits[-1]
              hFinal = ROOT.TVector3( hit.GetStart()[0], hit.GetStart()[1], hit.GetStart()[2] )
              t_TMSDeath[0] = hFinal.x()
              t_TMSDeath[1] = hFinal.y()
              t_TMSDeath[2] = hFinal.z()

            tout.Fill()

        ient += 1

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    #parser.add_option('--outfile', help='Output file name', default="out.root")
    #parser.add_option('--topdir', help='Input file top directory', default="")
    #parser.add_option('--nmax', help='Maximum number of events', default="99999999999999")
    #parser.add_option('--maxfiles', help='Maximum number of events', default="99999")
    #parser.add_option('--first_run', type=int, help='First run number', default=0)
    #parser.add_option('--last_run', type=int, help='Last run number', default=0)
    #parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    #parser.add_option('--grid', action='store_true', help='grid mode', default=False)
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

    t_lepTMS = array('f',3*[0.0])
    tout.Branch('lepTMS',t_lepTMS,'leTMS[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_TMSKE = array('f',[0])
    tout.Branch('TMSKE',t_TMSKE,'TMSKE/F')
    t_lepE = array('f',[0])
    tout.Branch('lepE',t_lepE,'lepE/F')

    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonExitKE = array('f',[0])
    tout.Branch('muonExitKE',t_muonExitKE,'muonExitKE/F')

    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muonStart = array('i',[0])
    tout.Branch('muonStart',t_muonStart,'muonStart/I')

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

    #nfiles = 0
    #for run in glob.glob(args.topdir+"/edep*.root"):
      #print nfiles
      #print maxfiles
      #if nfiles > maxfiles: 
        #break

      #fname = run
    fname = input_file
    print "Adding "+fname+" to TChain..."
    tgeo = None

    tf = ROOT.TFile.Open( fname )
    tf.MakeProject("EDepSimEvents","*","RECREATE++")

    # add it to the tchain
    events.Add( fname )
    dspt.Add( fname )

    #if tgeo is None: # first OK file, get geometry
    tgeo = tf.Get("EDepSimGeometry")
    tf.Close() # done with this one

    #nfiles += 1

    #print "Running on ",nfiles," root files..."
    loop( events, dspt, tgeo, tout )

    fout.cd()
    tout.Write()

