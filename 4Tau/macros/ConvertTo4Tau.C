void ConvertTo4Tau(double ma = 15,
		   double BRTo2Mu2Tau = 0.5e-3) {

  double m_tau = 1.777;
  double m_mu = 0.106;

  double denominator = m_mu*m_mu*TMath::Sqrt(1-4*m_mu*m_mu/(ma*ma));
  double numerator   = m_tau*m_tau*TMath::Sqrt(1-4*m_tau*m_tau/(ma*ma));

  double BRTo4Tau = BRTo2Mu2Tau*0.5*numerator/denominator;
  std::cout << std::endl;
  std::cout << "B(H->aa->2mu2tau) = " << BRTo2Mu2Tau << std::endl;
  std::cout << "B(H->aa->4tau)    = " << BRTo4Tau << std::endl;

  
}
