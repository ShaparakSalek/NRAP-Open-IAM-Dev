NRAP-Open-IAM
=============
This is the prototype NRAP Phase II Integrated Assessment Model (NRAP-Open-IAM).
To use the NRAP-Open-IAM, download and extract the folders, then open
the User_Guide.pdf in the root folder.  To download the NRAP-Open-IAM click on
the download icon in the lower right section of the header information
next to the Clone button.
The download icon looks like: ![gitlab download image](https://gitlab.com/NRAP/OpenIAM/-/wikis/uploads/ec2d7d5cb2b0c45d5b61e69acccfc8cc/download_icon.PNG)

From the menu that appears select the preferred compressed folder
data type. The download will have the current repository hash as part of the
folder name, feel free to rename the folder something simple
such as "OpenIAM".

Additionally, the archive of the tool can be obtained through
NETL's Energy Data eXchange website: https://edx.netl.doe.gov/dataset/nrap-open-source-iam
by requesting an access through e-mail addressed to NRAP@netl.doe.gov.

Docker
-------
We provide a Docker image that runs the GUI so that NRAP-Open-IAM can be tried out without installation.
Docker will need to be installed on your computer: https://www.docker.com/products/docker-desktop
Ensure that Docker has been started and you have signed into your Docker account.

Then, on MacOSX/Linux, run: 

xhost + 127.0.0.1

docker run -e DISPLAY=host.docker.internal:0 -v $(pwd)/output:/output dharp/openiam-gui:beta

Documentation
-------------
The NRAP-Open-IAM documentation can be accessed at https://nrap-iam.gitlab.io/UQ_example_setup

Feedback
---------
As this is a prototype of software being actively developed, we
are seeking any feedback or bug reports.  An online feedback form
can be found here: https://docs.google.com/forms/d/e/1FAIpQLSed5mcX0OBx1dLNmYGbmS4Vfc0mdOLapIzFqw-6vHoho9B19A/viewform?usp=sf_link

Feedback can also be emailed to the OpenIAM project at NRAP@netl.doe.gov, or the lead developer
at Veronika.Vasylkivska@netl.doe.gov or any other
member of the development team.

Questions can also be asked through the forum available at NETL's Energy Data eXchange website
https://edx.netl.doe.gov/organization/forum/nrap-tools.

Disclaimer
-----------
This software was prepared as part of work sponsored by an agency
of the United States Government. Neither the United States Government
nor any agency thereof, nor any of their employees, makes any warranty,
express or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information,
apparatus, product, or process disclosed, or represents that its use
would not infringe privately owned rights. Reference therein to any
specific commercial product, process, or service by trade name,
trademark, manufacturer, or otherwise does not necessarily constitute
or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof. The views and opinions of
authors expressed therein do not necessarily state or reflect those
of the United States Government or any agency thereof.
