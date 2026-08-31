import ROOT
ROOT.gROOT.SetBatch(True)

def make_cb(n, mean, sigma, alpha, norm):
    crystalball = ROOT.TF1("crystalball", "crystalball", min_x, max_x)

    crystalball.SetParameter(0, n) # n
    crystalball.SetParameter(1, mean)
    crystalball.SetParameter(2, sigma)
    crystalball.SetParameter(3, alpha) # alpha
    crystalball.SetParameter(4, norm) # norm
    
    return crystalball

min_x = -5
max_x = 2
gauss = ROOT.TF1("gauss", "gaus", min_x, max_x)
gauss.SetLineColor(ROOT.kRed)
gauss.SetLineStyle(2)
gauss.SetLineWidth(2)

norm = 1
mean = 0
sigma = 0.5
gauss.SetParameter(0, norm)
gauss.SetParameter(1, mean)
gauss.SetParameter(2, sigma)
#crystalball.SetParameters(2, 1, 1, 2, 3);

n = 1
crystalball = make_cb(n, mean, sigma, 0.5, 4)
crystalball.SetLineColor(ROOT.kBlack)
crystalball.SetLineWidth(2)
crystalball2 = make_cb(n, mean, sigma, 1, 1)
crystalball2.SetLineStyle(7)
crystalball2.SetLineColor(ROOT.kBlue)
crystalball2.SetLineWidth(2)

leg = ROOT.TLegend(0.15,0.55,0.60,0.85)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(gauss, "Gaussian, #bar{x}=0, #sigma=0.5", "L")
leg.AddEntry(crystalball, "n = 4, #alpha = 0.5", "L")
leg.AddEntry(crystalball2, "n = 1, #alpha = 1", "L")


canvas = ROOT.TCanvas()
crystalball.SetTitle("Crystal Ball Comparison")
crystalball.GetXaxis().SetTitle("X")
crystalball.GetYaxis().SetTitle("F(X)")

crystalball.Draw()
crystalball2.Draw("same")
gauss.Draw("same")
leg.Draw()
canvas.Print("crystal_ball.png")


