% Use matlab to run this file
% This code converts the 2D mesh file from matlab to julia format. To convert 
% the mesh, the mesh is first created in Gmsh and then exported in matlab format.
% Then, run this script to convert from '.m' extension to '.jl'.

clc; clear;

run('heat_conduction2D.m');
fid = fopen('heat_conduction2D.jl', 'w');

%% sort the msh.POS matrix in specific order as given below
[~, perm] = sortrows(round(msh.POS, 3), [3, 1, 2]);
msh.POS = msh.POS(perm, :);

%% build mapping from old to new indices
new_id = zeros(length(perm), 1);
for new_idx = 1:length(perm)
    old_idx = perm(new_idx);
    new_id(old_idx) = new_idx;
end

%% update QUADS connectivity (first 4 columns)
updated_QUADS = msh.QUADS;
for i = 1:size(msh.QUADS, 1)
    for j = 1:4
        updated_QUADS(i, j) = new_id(msh.QUADS(i, j));
    end
end
msh.QUADS = updated_QUADS;

%% initialize the msh dict
fprintf(fid, 'msh = Dict{Symbol, Any}()\n\n');

%% write the number of nodes in the mesh
fprintf(fid, 'msh[:nbNod] = %d\n\n', msh.nbNod);

%% write the maxmimum and minimum position in the mesh
fprintf(fid, 'msh[:max] = [%f, %f, %f]\n', msh.MAX);
fprintf(fid, 'msh[:min] = [%f, %f, %f]\n\n', msh.MIN);

%% write the POS array
fprintf(fid, 'msh[:POS] = [%0.12f %0.12f %0.12f\n', msh.POS(1,:));
for i=2:size(msh.POS,1)-1
    fprintf(fid, '\t\t\t  %0.12f %0.12f %0.12f\n', msh.POS(i,:));
end
fprintf(fid, '\t\t\t  %0.12f %0.12f %0.12f]\n\n', msh.POS(end,:));

%% write the QUADS array
fprintf(fid, 'msh[:QUADS] = [%d %d %d %d %d\n', msh.QUADS(1,1:end-1), msh.QUADS(1,end));
for i=2:size(msh.QUADS,1)-1
    fprintf(fid, '\t\t\t    %d %d %d %d %d\n', msh.QUADS(i,1:end-1), msh.QUADS(i,end));
end
fprintf(fid, '\t\t\t    %d %d %d %d %d]\n', msh.QUADS(end,1:end-1), msh.QUADS(end,end));

fclose(fid);