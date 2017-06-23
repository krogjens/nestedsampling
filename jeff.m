function [theta]=jeff(u,thmin,thmax)

theta=thmin*exp(u*log(thmax/thmin));