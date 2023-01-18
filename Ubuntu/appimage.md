# How to use AppImage in Ubuntu 22.04

1. Download AppImage file, e.g., `softwarename.AppImage`. 
1. Go to Download directry in terminal. 
    - `chmod a+x softwarename.AppImage`
    - Test launch with `./software.AppImage`
      - If core dump, try it with `--in-process-gpu` option. 

## Add to launcher
1. Move AppImage file to certain directry, e.g, `/opt/softwarename/`. 
1. Make `softwarename.desktop` file in `~/.local/share/applications/`. 
    - Below, the content of `mendeley-reference-manager.desktop` for example. 
    - Please put image file, e.g., `image.png`, in the path of Icon field if you need. 

`~/.local/share/applications/mendeley-reference-manager.desktop`
```
[Desktop Entry]
Name=Mendeley Referenced Manager
Comment=Mendeley Reference Manager is software for managing and sharing research papers
Exec=/opt/mendeley-reference-manager/mendeley-reference-manager-2.83.0-beta.0-x86_64.AppImage --in-process-gpu
Icon=/opt/mendeley-reference-manager/mendeley-reference-manager.png
Terminal=false
Type=Application
Categories=Education;Literature;Qt;
```

1. Check launcher and if there is the application of the AppImage. 

## Update AppImage
1. Download the latest AppImage file. 
1. Move AppImage file to certain directry, e.g, `/opt/softwarename/`. 
1. Update `softwarename.desktop`, i.e., rewrite the Exec path. 
