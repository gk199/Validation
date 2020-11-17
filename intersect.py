import numpy as np
import math
# this is a correction function for the eta phi value of the intersection of a particle not resulting from the IP       
# used since the b quark is matched to the TP based on eta phi, but the b quark results from a LLP decay so has a different vertex                
def intersect(vx, vy, vz, px, py, pz):
  lightSpeed = 29979245800 # cm/sec
  radius = 179 # 130cm for calorimeters (ECAL + HCAL)
  length = 388 # 300cm for calorimeters (ECAL + HCAL)
  p = math.sqrt(px*px + py*py + pz*pz)
  # First work out intersection with cylinder (barrel)        
  a = (px*px + py*py)*lightSpeed*lightSpeed/(p*p)
  b = 2*(vx*px + vy*py)*lightSpeed/p
  c = (vx*vx + vy*vy) - radius*radius
  sqrt_disc = math.sqrt(b*b - 4*a*c)
  tCircle1 = (-b + sqrt_disc)/(2*a)
  tCircle2 = (-b - sqrt_disc)/(2*a)
  # If intersection in the past it doesn't count         
  if (tCircle1 < 0):
     tCircle1 = 1e9
  if (tCircle2 < 0):
     tCircle2 = 1e9
  # If the intsersection occurs outside the barrel length it doesn't count                       
  zPosCircle1 = tCircle1*(pz/p)*lightSpeed + vz
  zPosCircle2 = tCircle2*(pz/p)*lightSpeed + vz
  if (zPosCircle1 > length):
     tCircle1 = 1e9;
  if (zPosCircle2 > length):
     tCircle2 = 1e9;
  # Now work out if it intersects the endcap                      
  tPlane1 = (length-vz)*p/(pz*lightSpeed)
  tPlane2 = (-length-vz)*p/(pz*lightSpeed)
  # If intersection in the past it doesn't count                     
  if (tPlane1 < 0):
     tPlane1 = 1e9
  if (tPlane2 < 0):
     tPlane2 = 1e9
  xPosPlane1 = tPlane1*(px/p)*lightSpeed + vx
  yPosPlane1 = tPlane1*(py/p)*lightSpeed + vy
  xPosPlane2 = tPlane2*(px/p)*lightSpeed + vx
  yPosPlane2 = tPlane2*(py/p)*lightSpeed + vy
  # If the intsersection occurs outside the endcap radius it doesn't count     
  if (math.sqrt(xPosPlane1*xPosPlane1 + yPosPlane1*yPosPlane1) > radius):
     tPlane1 = 1e9
  if (math.sqrt(xPosPlane2*xPosPlane2+yPosPlane2*yPosPlane2) > radius):
     tPlane2 = 1e9
  # Find the first intersection                          
  tInter = min(tCircle1,tCircle2,tPlane1,tPlane2)
  # Return 1000,1000 if not intersection with barrel or endcap             
  #if (tInter > 1E6)
  #  {
  #    etaphi.push_back(1000);
  #    etaphi.push_back(1000);
  #    return etaphi;
  #  }
  # Find position of intersection                          
  xPos = tInter*(px/p)*lightSpeed + vx
  yPos = tInter*(py/p)*lightSpeed + vy
  zPos = tInter*(pz/p)*lightSpeed + vz
  # Find eta/phi of intersection                          
  phi = math.atan2(yPos,xPos) # return the arc tan in radians
  theta = math.acos(zPos/math.sqrt(xPos*xPos + yPos*yPos + zPos*zPos))
  eta = -math.log(math.tan(theta/2.))
  return eta,phi

def deltaPhi(phi1, phi2):
  # Delta phi calculation, using 0 to 2pi
  if (phi1<0):
      phi1+=2.*math.pi
  if (phi2<0):
      phi2+=2.*math.pi
  result = phi1 - phi2
  #if(fabs(result) > 9999) 
  #  return result;
  while (result > math.pi):
    result -= 2.*math.pi
  while (result <= -math.pi):
    result += 2.*math.pi
  return result

def deltaR(eta1, phi1, eta2, phi2):
  deta = eta1 - eta2
  dphi = deltaPhi(phi1, phi2)
  return math.sqrt(deta*deta + dphi*dphi)

#     // find partons - quarks (d, y, s, c, b) or gluon from the gen particles. From QCD or LLP decay. top quark is unstable and heavier
#     // add detector cuts as well - done here before partons are saved to Lorentz vector (require high enough pt). Eta cut is also done
#     if ( (abs(pdgId)>=1 && abs(pdgId)<=5) || abs(pdgId)==21 ) {
#       partons.push_back(p4);
#       //std::cout << "pdg ID: " << it.pdgId() << std::endl;
#       //std::cout << "particle vertex: " << it.vertex() << std::endl;
#       //std::cout << "r position: " << it.vertex().R() << std::endl;
#       //std::cout << "x, y, z position: " << it.vertex().X() << it.vertex().Y() << it.vertex().Z() << std::endl;
#       double px = it.px();
#       double py = it.py();
#       double pz = it.pz();
#       double vx = it.vertex().X();
#       double vy = it.vertex().Y();
#       double vz = it.vertex().Z();
#       //       std::cout << intersect(vx, vy, vz, px, py, pz)[0] << "  eta, phi of parton " << intersect(vx, vy, vz, px, py, pz)[1] << std::endl;
#       // now add detector cuts of pt and eta, before the partons eta and phi vectors are saved. These are vectors used for delta R matching
#       if (it.pt()>20. && fabs(intersect(vx, vy, vz, px, py, pz)[0])<2.5 )
#	 {
#	   partonseta.push_back(intersect(vx, vy, vz, px, py, pz)[0]);
#	   partonsphi.push_back(intersect(vx, vy, vz, px, py, pz)[1]);
#	 }
#       // check that the LLP decays in the HCAL by requiring quark gen vertex to be in the HCAL volume, and then save these events to the Lorentz vector partonsHCAL
#       // HE region 3.88 - 5.68 m and out to 2.95 m radius
#       // HB region z below 3.88m, radius 1.79 - 2.95m
#       double radius = 0.;
#       radius = sqrt(it.vertex().X()*it.vertex().X()+it.vertex().Y()*it.vertex().Y());
#       // vertex x, y, z are in cm, radius is in cm, define the detector regions by cm as well. HCAL barrel 179-295 cm, HE 388-568
#       if ( abs(it.vertex().Z()) < 568 && radius < 295 && ( abs(it.vertex().Z()) > 388 || radius > 179 ) ) {
#	 partonsHCAL.push_back(p4);
#	 if (it.pt()>20. && fabs(intersect(vx, vy, vz, px, py, pz)[0])<2.5 )
#	   {
#	     partonsetaHCAL.push_back(intersect(vx, vy, vz, px, py, pz)[0]);
#	     partonsphiHCAL.push_back(intersect(vx, vy, vz, px, py, pz)[1]);
#	   }
#	 //	 std::cout << "radius: " << radius << std::endl;
#	 //	 std::cout << "particle partons vertex: " << it.vertex() << std::endl;
#	 //      std::cout << intersect(vx, vy, vz, px, py, pz)[0] << "  eta, phi of parton " << intersect(vx, vy, vz, px, py, pz)[1] << std::endl;
#       }
#     }
