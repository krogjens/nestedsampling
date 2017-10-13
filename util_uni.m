function [theta]=uni(u,theta_min,theta_max)

theta=(theta_max - theta_min)*u +theta_min;
