
Double_t rms90(TH1 *h) {

  TAxis *axis = h->GetXaxis();
  Int_t nbins = axis->GetNbins();
  Int_t imean = axis->FindBin(h->GetMean());
  //  Double_t entries = 0.9 * h->GetEntries();
  Double_t entries = 0.9 * h->GetSumOfWeights();
  Double_t w = h->GetBinContent(imean);
  Double_t x = h->GetBinCenter(imean);

  Double_t sumw = w;
  Double_t sumwx = w*x;
  Double_t sumwx2 = w*x*x;

  for (Int_t i=1;i<nbins;i++) {
    if ( imean-i > 0) {
      w = h->GetBinContent(imean-i);
      x = h->GetBinCenter(imean-i);
      sumw += w;
      sumwx += w*x;
      sumwx2 += w*x*x;
    }
    if ( imean+i <= nbins) {
      w = h->GetBinContent(imean+i);
      x = h->GetBinCenter(imean+i);
      sumw += w;
      sumwx += w*x;
      sumwx2 += w*x*x;
    }
    //    printf(" i=%d: sumw =%g - entries =%g \n" , i, sumw, entries ) ;
    if (sumw > entries) {
      //      printf(" ---- breaking : sumw =%g - entries =%g,  i=%d" , sumw, entries, i ) ;
      break;
    }
  }
  x = sumwx/sumw;
  Double_t rms2 = TMath::Abs(sumwx2/sumw -x*x);
  Double_t result = TMath::Sqrt(rms2);

//   printf("RMS of central 90%% = %g, RMS total = %g\n, Mean90 = %g , Mean total = %g \n " ,
	//  result, h->GetRMS() , x , h->GetMean() );

  return result ;
}

double DoubleSidedCrystalballFunction(double *x, double *par)
{
  double alpha_l = par[0];
  double alpha_h = par[1];
  double n_l     = par[2];
  double n_h     = par[3];
  double mean	= par[4];
  double sigma	= par[5];
  double N	= par[6];
  double t = (x[0]-mean)/sigma;
  double result;
  double fact1TLessMinosAlphaL = alpha_l/n_l;
  double fact2TLessMinosAlphaL = (n_l/alpha_l) - alpha_l -t;
   double fact1THihgerAlphaH = alpha_h/n_h;
   double fact2THigherAlphaH = (n_h/alpha_h) - alpha_h +t;

   if (-alpha_l <= t && alpha_h >= t)
     {
       result = exp(-0.5*t*t);
     }
   else if (t < -alpha_l)
     {

       result = exp(-0.5*alpha_l*alpha_l)*pow(fact1TLessMinosAlphaL*fact2TLessMinosAlphaL, -n_l);

     }
   else if (t > alpha_h)
     {
       result = exp(-0.5*alpha_h*alpha_h)*pow(fact1THihgerAlphaH*fact2THigherAlphaH, -n_h);

     }
   return N*result;
}
