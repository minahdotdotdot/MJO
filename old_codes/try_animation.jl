using PyPlot
using PyCall
using Base64
@pyimport matplotlib.animation as anim

#Construct Figure and Plot Data
fig = figure("MyFigure",figsize=(5,5))
ax = PyPlot.axes(xlim = (0,10),ylim=(0,10))
global line1 = ax[:plot]([],[],"r-")[1]
global line2 = ax[:plot]([],[],"g-")[1]
global line3 = ax[:plot]([],[],"b-")[1]

# Define the init function, which draws the first frame (empty, in this case)
function init()
    global line1
    global line2
    global line3
    line1[:set_data]([],[])
    line2[:set_data]([],[])
    line3[:set_data]([],[])
    return (line1,line2,line3,Union{})  # Union{} is the new word for None
end

# Animate draws the i-th frame, where i starts at i=0 as in Python.
function animate(i)
    global line1
    global line2
    global line3
    x = (0:i)/10.0
    line1[:set_data](x,x)
    line2[:set_data](1 .+x,x)
    line3[:set_data](2 .+x,x)
    return (line1,line2,line3,Union{})
end

# Create the animation object by calling the Python function FuncAnimaton
myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20);

# Convert it to an MP4 movie file and saved on disk in this format.
myanim[:save]("3Lines.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"]);
#myanim[:save]("test1.mp4", bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
    
# Function for creating an embedded video given a filename
function html_video(filename)
    open(filename) do f
        base64_video = base64encode(f)
        """<video controls src="data:video/x-m4v;base64,$base64_video">"""
    end
end

# Display the movie in a Julia cell as follows. Note it has animation controls for the user.
display("text/html", html_video("3Lines.mp4"))
#display("text/html", html_video("test1.mp4"))




