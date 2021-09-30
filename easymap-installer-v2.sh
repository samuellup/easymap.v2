#!/bin/bash
# This is the Easymap ato-installation script, its meant for facilitating the installation of Easymap for the general user in different Unix distributions. 
# Administraror privileges are required for the execution of this installer. 

# Links to the main Easymap repositories
git_address="https://github.com/samuellup/easymap.v2.git"     
sf_download="wget -O easymap-v2.zip https://sourceforge.net/projects/easymap-v2/files/latest/download"
# Default access port: 8100, can be changed adding a port number as an argument when running the script
if ! [ $1 ]; then
	port=8100

elif [ "$1" -ge 8100 ] && [ "$1" -le 8200 ]; then
	port=$1

else
	echo 'Please choose a port number between 8100 and 8200.'
	exit
fi

# First the user selects the Operative System used with a  simple menu. Then the script installs different dependencies in each case and clones 
# easymap from the main repository. Finally, it runs the install.sh file inside the Easymap folder for a default installation using port 8100
# as the access point to the web interface. 

echo '
Welcome to the Easymap installer, please select your Operating System:

'
PS3='

Enter a number: '

distrs='Ubuntu_20 Ubuntu_18 Ubuntu_16 Ubuntu_14 Linux_AMI Linux_Redhat Quit'
select dis in $distrs
do
	if [ $dis == 'Quit' ]
	then
		echo "
Installation cancelled. 
			" 
		break
	else
		echo '
Beggining Easymap installation in' $dis'. Please wait for the process to finish, it could take up to 30 minutes.
'
		sleep 4s
	fi

	if [ $dis == 'Ubuntu_20' ]
	then
		apt-get update
		apt-get install build-essential zlib1g-dev libbz2-dev git wget tar zip liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev make -y
		if [ -d easymap.v2 ]; then rm -rf easymap.v2; fi
		#git clone $git_address
		$sf_download
		unzip easymap-v2.zip
		chmod -R 755 easymap.v2
		cd easymap.v2
		./install.sh server $port
	fi

	if [ $dis == 'Ubuntu_18' ]
	then
		apt-get update  
		apt-get install build-essential zlib1g-dev libbz2-dev git wget tar zip liblzma-dev libncurses5-dev libncursesw5-dev libssl1.0-dev -y
		if [ -d easymap.v2 ]; then rm -rf easymap.v2; fi
		#git clone $git_address
		$sf_download
		unzip easymap-v2.zip
		chmod -R 755 easymap.v2
		cd easymap.v2
		./install.sh server $port
	fi

    if [ $dis == 'Ubuntu_16' ] || [ $dis == 'Ubuntu_14' ]
    then
		apt-get update
		apt-get install build-essential zlib1g-dev libbz2-dev git wget tar zip liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev -y
		if [ -d easymap.v2 ]; then rm -rf easymap.v2; fi
		#git clone $git_address
		$sf_download
		unzip easymap-v2.zip
		chmod -R 755 easymap.v2
		cd easymap.v2
		./install.sh server $port
    fi

    if [ $dis == 'Linux_AMI' ]
    then
        yum groupinstall -y "Development Tools"
		yum groupinstall -y "Development Libraries"
		yum install -y wget zlib-devel bzip2-devel git ncurses-devel ncurses openssl-devel
		yum install -y xz-devel
		if [ -d easymap.v2 ]; then rm -rf easymap.v2; fi
		#git clone $git_address
		$sf_download
		unzip easymap-v2.zip
		chmod -R 755 easymap.v2
		cd easymap.v2
		./install.sh cli
		sudo -H -u $SUDO_USER bash -c "nohup ./src/Python-2.7.18/.localpython/bin/python -m CGIHTTPServer $port" &
	fi

    if [ $dis == 'Linux_Redhat' ]
    then
		yum groupinstall -y "Development Tools"
		yum groupinstall -y "Development Libraries"
		yum install -y wget zlib-devel bzip2-devel git ncurses-devel ncurses openssl-devel
		yum install -y xz-devel
		yum install -y curl-devel
		yum install -y curl
		if [ -d easymap.v2 ]; then rm -rf easymap.v2; fi
		#git clone $git_address
		$sf_download
		unzip easymap-v2.zip
		chmod -R 755 easymap.v2
		cd easymap.v2
		./install.sh cli
		sudo -H -u $SUDO_USER bash -c "nohup ./src/Python-2.7.18/.localpython/bin/python -m CGIHTTPServer $port" &
    fi
	break
done

