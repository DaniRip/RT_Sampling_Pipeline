function [clusters, centroids] = kMeansClustering(mat_file,numClusters,numIterations)
    input = load(mat_file);
    input = input.guiInput;
    num_voxels = size(input.Dij,1);
    num_beamlets = size(input.Dij,2);
    target_dose = input.targetDose;
    beam_width = input.beamWidth;
    d_beam_indicies = input.beamIndicies;

    avgPoints = rand(numClusters,num_beamlets);
    for j = 1:num_beamlets
        avgPoints(:,j) = avgPoints(:,j)*(max(input.Dij(:,j))-min(input.Dij(:,j)))+min(input.Dij(:,j));
    end
    
    for i = 1:numClusters
        j = ceil(rand*num_voxels);
        while sum(ismember(avgPoints,input.Dij(j,:),'rows')) ~= 0
            j = ceil(rand*num_voxels);
        end
        avgPoints(i,:) = input.Dij(j,:);
    end
    
    for iter = 1:numIterations
%         disp(strcat('Iteration:',{' '},string(iter)));
        dataSetAssignments = [input.Dij ones(num_voxels,1)];
        for i = 1:size(dataSetAssignments,1)
           minDist = norm(dataSetAssignments(i,1:num_beamlets).' - avgPoints(1,:).');
           minJ = 1;
           for j = 1:size(avgPoints,1)
               dist = norm(dataSetAssignments(i,1:num_beamlets).' - avgPoints(j,:).');
               if dist <= minDist
                   minJ = j;
                   minDist = dist;
               end
           end
           dataSetAssignments(i,num_beamlets+1) = minJ;
        end
        for i = 1:numClusters
            splitSet(:,:,i) = {dataSetAssignments(dataSetAssignments(:,num_beamlets+1)==i,1:num_beamlets)};
        end
        
        
        for i = 1:numClusters
            avg = mean(splitSet{i},1);
            avgPoints(i,:) = avg(1:num_beamlets);
        end
    end
    
    clusters = splitSet;
    centroids = avgPoints;
    
end