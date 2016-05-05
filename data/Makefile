url = https://s3.amazonaws.com/3kricegenome/reduced

NB-core_v4 : 
	wget -nv -N ${url}/$@.map.gz
	wget -nv -N ${url}/$@.ped.gz
	sha1sum -c $@.shasum


3krg_filt_snp_v4 : 
	wget -nv -N ${url}/$@.map.gz
	wget -nv -N ${url}/$@.ped.gz
	sha1sum -c $@.shasum
