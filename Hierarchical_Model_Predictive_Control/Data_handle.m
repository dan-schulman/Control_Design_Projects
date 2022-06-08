classdef Data_handle
    methods(Static)
        function [data] = dataCollect(filename)
            % Collecting Data
            delimiterIn = '\t';
            headerlinesIn = 14;
            A = importdata(filename,delimiterIn,headerlinesIn);

            rawTimeStamp = A.data(:,1);
            rawData = A.data(:,2);
            timeStamp = A.data(:,1)-A.data(1,1);
            %Normalizing data
            normTimeStamp = A.data(:,1)-A.data(1,1);
            normData = A.data(:,2)-mean(A.data(1:5,2));

            %Data Filtering
            filteredData = smoothdata(rawData,'gaussian',80);

            %Saving data to struct which is returned
            data.timeStamp = timeStamp;
            data.rawData = rawData;
            data.normTimeStamp = normTimeStamp;
            data.normData = normData;
            data.filteredData = filteredData;
        end
        
        function [h,v] = readtsv(dataDir)
            currentDir = pwd;
            myFiles = dir(fullfile(dataDir,'*.tsv'));
            fileNames = string({myFiles.name});
            cd(dataDir);
            % Head Data
            h.ax = Data_handle.dataCollect(fileNames(contains(fileNames,'C0008')==true));
            h.ay = Data_handle.dataCollect(fileNames(contains(fileNames,'C0009')==true));
            h.az = Data_handle.dataCollect(fileNames(contains(fileNames,'C0010')==true));
            h.wx = Data_handle.dataCollect(fileNames(contains(fileNames,'C0011')==true));
            h.wy = Data_handle.dataCollect(fileNames(contains(fileNames,'C0012')==true));
            h.wz = Data_handle.dataCollect(fileNames(contains(fileNames,'C0013')==true));

            %Vehicle Data
            v.ax = Data_handle.dataCollect(fileNames(contains(fileNames,'C0028')==true));
            v.ay = Data_handle.dataCollect(fileNames(contains(fileNames,'C0029')==true));
            v.az = Data_handle.dataCollect(fileNames(contains(fileNames,'C0030')==true));
            v.wx = Data_handle.dataCollect(fileNames(contains(fileNames,'C0031')==true));
            v.wy = Data_handle.dataCollect(fileNames(contains(fileNames,'C0032')==true));
            v.wz = Data_handle.dataCollect(fileNames(contains(fileNames,'C0033')==true));
            cd(currentDir)
        end
        function [ax,ay,az,wx,wy,wz,timeStop] = setAcc(h)
            ax = timeseries(h.ax.filteredData, h.ax.timeStamp);
            ay = timeseries(h.ay.filteredData, h.ay.timeStamp);
            az = timeseries(h.az.filteredData, h.az.timeStamp);

            wx = timeseries(h.wx.filteredData, h.wx.timeStamp);
            wy = timeseries(h.wy.filteredData, h.wy.timeStamp);
            wz = timeseries(h.wz.filteredData, h.wz.timeStamp);
            timeStop = h.ax.timeStamp(length(h.ax.timeStamp)-1);
        end
        
        function [ax,ay,az,wx,wy,wz,timeStop] = setVehAcc(v)
            ax = timeseries(v.ax.filteredData, v.ax.timeStamp);
            ay = timeseries(v.ay.filteredData, v.ay.timeStamp);
            az = timeseries(v.az.filteredData, v.az.timeStamp);

            wx = timeseries(v.wx.filteredData, v.wx.timeStamp);
            wy = timeseries(v.wy.filteredData, v.wy.timeStamp);
            wz = timeseries(v.wz.filteredData, v.wz.timeStamp);
            timeStop = v.ax.timeStamp(length(v.ax.timeStamp)-1);
        end
        
        function [ax_vis,ay_vis,az_vis,wx_vis,wy_vis,wz_vis] = setVisualAcc(conflict,ax,ay,az,wx,wy,wz)
            switch conflict
                case 1
                    ax_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                    ay_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                    az_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                    wx_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                    wy_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                    wz_vis = timeseries(zeros(1,length(ax.Time)),zeros(1,length(ax.Time)));
                case 0
                    ax_vis = ax;
                    ay_vis = ay;
                    az_vis = az;
                    wx_vis = wx;
                    wy_vis = wy;
                    wz_vis = wz;
            end
        end
        


        
    end
end