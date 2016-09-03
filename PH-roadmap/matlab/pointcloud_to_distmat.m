function [ D] = pointcloud_to_distmat( filename )
%Input: text file with point cloud, one point per line
%Output: text file with distance matrix, one row per line

%Nina Otter, Oxford 2015

A=load(filename);
size(A,1);
size(A,2);
dist=pdist(A);

D=squareform(dist);
end


