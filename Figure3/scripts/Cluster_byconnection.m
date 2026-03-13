function [node_groupid,totalgnum,changedleaforder_bygroups] = Cluster_byconnection(distm,distthresh,optleaforder)
% distm: distance matrix
% distthresh: threshold
% optleaforder: optimal leaf order 
% seperate the connections by distance therehold
datanum=size(distm,1);

if nargin<3
    optleaforder=1:datanum;
end
connect_m=distm(optleaforder,optleaforder)<distthresh; % find connection matrix by threshold (binaryzation)

%figureimagesc(connect_m)

totalgnum=1; % total nember of group
node_used=zeros(1,datanum); % already used for search
node_groupid_aloneleaf=zeros(1,datanum); % record group id of each node

currentnode=1; 
currentgroupmembers=currentnode;

while min(node_groupid_aloneleaf)==0 % not all of them are assigned
    newmembers=find(connect_m(currentnode,:)==1 & node_groupid_aloneleaf==0); % find members connect with current node
    node_groupid_aloneleaf(newmembers)=totalgnum; % set group number of new member as the same group with current node
    node_used(currentnode)=1; % current node is used 
    currentgroupmembers=[currentgroupmembers,newmembers]; % add newmember to current group
    nexttouse=min(currentgroupmembers(node_used(currentgroupmembers)==0)); % 下一个node，在现在group中且不在node_used中. in next iteration, find node connect with next node

    if ~isempty(nexttouse) % 如果存在next node，设为current node
        currentnode=nexttouse;
    else
        totalgnum=totalgnum+1;
        if min(node_groupid_aloneleaf)~=0
            totalgnum = totalgnum-1;
        end
        currentnode=min(find(node_groupid_aloneleaf==0));
        % switch to another group
        currentgroupmembers=currentnode;
        % add current node to current group
    end
end



node_groupid=0*node_groupid_aloneleaf;
node_groupid(optleaforder)=node_groupid_aloneleaf; % group id, reorder with optimal leaf order



changedleaforder_bygroups=[];
for i=1:totalgnum
    membersinorder=find(node_groupid(optleaforder)==i); % find node in each group
    changedleaforder_bygroups=[changedleaforder_bygroups,optleaforder(membersinorder)]; % 将这些cluster排列成order
end

end



